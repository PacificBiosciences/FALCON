"""Exact copy of falcon_unzip/tasks/snakemake.py
TODO: Consolidate.
"""
from __future__ import absolute_import


from future.utils import viewitems
from future.utils import itervalues

from builtins import object
import json
import os
import re


def find_wildcards(pattern):
    """
    >>> find_wildcards('{foo}/{bar}')
    ['bar', 'foo']
    """
    re_wildcard = re.compile(r'\{(\w+)\}')
    found = [mo.group(1) for mo in re_wildcard.finditer(pattern)]
    return list(sorted(found))

class SnakemakeRuleWriter(object):
    def legalize(self, rule_name):
        return self.re_bad_char.sub('_', rule_name, count=0)
    def unique_rule_name(self, basename):
        rule_name = basename
        if rule_name in self.rule_names:
            i = 1
            while rule_name in self.rule_names:
                rule_name = basename + str(i)
                i += 1
        self.rule_names.add(rule_name)
        return rule_name
    def write_dynamic_rules(self, rule_name, input_json, inputs, shell_template,
            parameters, wildcard_outputs, output_json):
        """Lots of conventions.
        input_json: should have a key 'mapped_inputs', which is a map of key->filename
          Those filenames will be symlinked here, according to the patterns in wildcard_inputs.
        shell_template: for the parallel task
        output_json: This will contain only key->filename, based on wildcard_outputs.
        inputs: These include non-wildcards too.
        (For now, we assume inputs/outputs is just one per parallel task.)
        """
        # snakemake does not like paths starting with './'; they can lead to mismatches.
        # So we run normpath everywhere.
        input_json = os.path.normpath(input_json)
        output_json = os.path.normpath(output_json)

        # snakemake cannot use already-generated files as dynamic outputs (the wildcard_inputs for the parallel task),
        # so we rename them and plan to symlink.
        wildcard_inputs = dict(inputs)
        nonwildcard_inputs = dict()
        for (key, fn) in list(viewitems(wildcard_inputs)):
            if '{' not in fn:
                del wildcard_inputs[key]
                nonwildcard_inputs[key] = fn
                continue
            dn, bn = os.path.split(wildcard_inputs[key])
            wildcard_inputs[key] = os.path.join(dn + '.symlink', bn)
        rule_name = self.unique_rule_name(rule_name)
        dynamic_output_kvs = ', '.join("%s=dynamic('%s')"%(k, os.path.normpath(v)) for (k, v) in viewitems(wildcard_inputs))
        dynamic_input_kvs =  ', '.join("%s=ancient(dynamic('%s'))"%(k, os.path.normpath(v)) for (k, v) in viewitems(wildcard_outputs))
        rule_parameters = {k: v for (k, v) in viewitems(parameters) if not k.startswith('_')}
        params = ','.join('\n        %s="%s"'%(k,v) for (k, v) in viewitems(rule_parameters))
        pattern_kv_list = list()
        for (name, wi) in viewitems(wildcard_inputs):
            fn_pattern = wi
            fn_pattern = fn_pattern.replace('{', '{{')
            fn_pattern = fn_pattern.replace('}', '}}')
            pattern_kv_list.append('%s="%s"' %(name, fn_pattern))
        wi_pattern_kvs = ' '.join(pattern_kv_list)

        rule = """
rule dynamic_%(rule_name)s_split:
    input:  %(input_json)r
    output: %(dynamic_output_kvs)s
    shell: 'python -m falcon_kit.mains.copy_mapped --special-split={input} %(wi_pattern_kvs)s'
"""%(locals())
        self.write(rule)

        input_wildcards = set() # Not sure yet whether input must match output wildcards.
        for wi_fn in itervalues(wildcard_inputs):
            found = find_wildcards(wi_fn)
            input_wildcards.update(found)
        wildcards = list(sorted(input_wildcards))
        params_plus_wildcards = {k: '{%s}'%k for k in wildcards}
        params_plus_wildcards.update(parameters)
        # The parallel script uses all inputs, not just wildcards.
        all_inputs = dict(wildcard_inputs)
        all_inputs.update(nonwildcard_inputs)
        self.write_script_rule(all_inputs, wildcard_outputs, params_plus_wildcards, shell_template, rule_name=None)

        wo_str_lists_list = ['%s=[str(i) for i in input.%s]' %(name, name) for name in list(wildcard_outputs.keys())]
        wo_pattern_kv_list = ['%s="%s"' %(name, os.path.normpath(patt)) for (name, patt) in viewitems(wildcard_outputs)]
        wo_str_lists_kvs = ',\n              '.join(wo_str_lists_list)
        wo_pattern_kvs =   ',\n              '.join(wo_pattern_kv_list)

        wildcards = list()
        for wi_fn in itervalues(wildcard_outputs):
            found = find_wildcards(wi_fn)
            if wildcards:
                assert wildcards == found, 'snakemake requires all outputs (and inputs?) to have the same wildcards'
            else:
                wildcards = found
        wildcards_comma_sep = ', '.join('"%s"' %k for k in wildcards)

        rule = '''
rule dynamic_%(rule_name)s_merge:
    input:  %(dynamic_input_kvs)s
    output: %(output_json)r
    run:
        snake_merge_multi_dynamic(output[0],
            dict(
              %(wo_str_lists_kvs)s
            ),
            dict(
              %(wo_pattern_kvs)s
            ),
            [%(wildcards_comma_sep)s] # all wildcards
        )
'''%(locals())
        self.write(rule)
    def write_script_rule(self, inputs, outputs, parameters, shell_template, rule_name):
        assert '_bash_' not in parameters
        first_output_name, first_output_fn = outputs.items()[0] # for rundir, since we cannot sub wildcards in shell
        if not rule_name:
            rule_name = os.path.dirname(first_output_fn)
        rule_name = self.unique_rule_name(self.legalize(rule_name))
        wildcard_rundir = os.path.normpath(os.path.dirname(first_output_fn)) # unsubstituted
        # We use snake_string_path b/c normpath drops leading ./, but we do NOT want abspath.
        input_kvs = ', '.join('%s=%s'%(k, snake_string_path(v)) for k,v in
                sorted(viewitems(inputs)))
        output_kvs = ', '.join('%s=%s'%(k, snake_string_path(v)) for k,v in
                sorted(viewitems(outputs)))
        rule_parameters = {k: v for (k, v) in viewitems(parameters) if not k.startswith('_')}
        #rule_parameters['reltopdir'] = os.path.relpath('.', wildcard_rundir) # in case we need this later
        params = ','.join('\n        %s="%s"'%(k,v) for (k, v) in viewitems(rule_parameters))
        shell = snake_shell(shell_template, wildcard_rundir)
        # cd $(dirname '{output.%(first_output_name)s}')
        rule = """
rule static_%(rule_name)s:
    input:  %(input_kvs)s
    output: %(output_kvs)s
    params:%(params)s
    shell:
        '''
outdir=$(dirname {output[0]})
#mkdir -p ${{outdir}}
cd ${{outdir}}
date

%(shell)s
date
'''
"""%(locals())
        self.write(rule)
    def __call__(self, inputs, outputs, parameters, shell_template, rule_name=None):
        self.write_script_rule(inputs, outputs, parameters, shell_template, rule_name)
    def __init__(self, writer):
        self.write = writer.write
        self.rule_names = set() # to ensure uniqueness
        self.re_bad_char = re.compile(r'\W')
        self.write("""
# THIS IS CURRENTLY BROKEN.
import json
import os
#import snakemake.utils

def snake_merge_dynamic_dict(reldir, input_fns, pattern, wildcards):
        '''Assume each wildcard appears at most once in the pattern.
        '''
        for k in wildcards:
            pattern = pattern.replace('{%s}' %k, '(?P<%s>\w+)' %k)
        re_dynamic = re.compile(pattern)
        mapped = list()
        for fn in input_fns:
            mo = re_dynamic.search(fn)
            assert mo, '{!r} did not match {!r}'.format(fn, re_dynamic.pattern)
            file_description = dict()
            file_description['wildcards'] = dict(mo.groupdict())
            file_description['fn'] = os.path.relpath(fn, reldir)
            mapped.append(file_description)
        return mapped

def snake_merge_multi_dynamic(output_fn, dict_of_input_fns, dict_of_patterns, wildcards):
        outdir = os.path.normpath(os.path.dirname(output_fn))
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        assert list(sorted(dict_of_input_fns.keys())) == list(sorted(dict_of_patterns.keys()))
        all_mapped = dict()
        for i in dict_of_patterns.keys():
            input_fns = dict_of_input_fns[i]
            pattern = dict_of_patterns[i]
            mapped = snake_merge_dynamic_dict(outdir, input_fns, pattern, wildcards)
            all_mapped[i] = mapped
        all_grouped = dict()
        for i, mapped in all_mapped.items():
            #print(i, mapped)
            for file_description in mapped:
                #print(file_description)
                #print(file_description['wildcards'])
                #print(list(sorted(file_description['wildcards'].items())))
                wildkey = ','.join('{}={}'.format(k,v) for k,v in sorted(file_description['wildcards'].items()))
                if wildkey not in all_grouped:
                    new_group = dict(
                        wildcards=dict(file_description['wildcards']),
                        fns=dict(),
                    )
                    all_grouped[wildkey] = new_group
                group = all_grouped[wildkey]
                wildcards = file_description['wildcards']
                assert wildcards == group['wildcards'], '{!r} should match {!r} by snakemake convention'.format(
                    wildcards, group['wildcards'])
                fn = file_description['fn']
                group['fns'][i] = fn
        ser = json.dumps(all_grouped, indent=2, separators=(',', ': ')) + '\\n'
        with open(output_fn, 'w') as out:
            out.write(ser)
""")
        prefix = """
shell.prefix('''
# Add -e vs. in falcon_unzip.
set -vex
hostname
pwd
''')
"""
        self.write(prefix)
class SnakemakeDynamic(object):
    """Not currently used."""
    def __init__(self, path):
        self.path = path
def snake_string_path(p):
    """normpath drops leading ./
    """
    if isinstance(p, SnakemakeDynamic):
        return "dynamic('{}')".format(
                os.path.normpath(p.path))
    else:
        return "'{}'".format(
                os.path.normpath(p))
def snake_shell(template, rundir):
    reltopdir = os.path.relpath('.', rundir)
    def makerel(mo):
        return os.path.join(reltopdir, mo.group(0))
    re_inout = re.compile(r'{(?:input|output)')
    return re_inout.sub(makerel, template, count =0)
