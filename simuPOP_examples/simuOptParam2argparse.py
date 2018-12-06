#!/usr/bin/env python

#
#
# To perform conversion:
#
# python simuOptParam2argparse.py filename start_lineno ending_lineno
#

import argparse


def convert_options(filename, startline, endline, name):
    with open(filename) as script:
        content = ''.join(script.readlines()[startline: endline])
    prefix = '''
class FakeSimuOptClass:
    def __init__(self):
        for name in [
            'valueNot',
            'valueOr',
            'valueAnd',
            'valueOneOf',
            'valueTrueFalse',
            'valueBetween',
            'valueGT',
            'valueGE',
            'valueLT',
            'valueLE',
            'valueEqual',
            'valueNotEqual',
            'valueIsInteger',
            'valueIsNum',
            'valueIsList',
            'valueValidDir',
            'valueValidFile',
            'valueListOf',
            'valueSumTo']:
            setattr(self, name, lambda *args, **kwargs: None)
simuOpt = FakeSimuOptClass()
'''
    # replace list
    replaces = [
    ]

    for old, new in replaces:
        content = content.replace(old, new)
    value = {}
    exec(prefix + content, value)

    print('import argparse')
    print('args = argparse.ArgumentParser()')
    for opt in value[name]:
        if 'name' not in opt:
            continue
        name = opt['name']
        default = "\n    default={!r},".format(opt['default']) if 'default' in opt else None
        type = opt['type'] if 'type' in opt else None
        help = opt['label'] if 'label' in opt else ''
        help = (help + '\n        ' + opt['description']) if 'description' in opt else help
        if 'type' in opt:
            if opt['type'] == 'filename':
                topt = ''
            elif not isinstance(opt['type'], (list, tuple)):
                if hasattr(opt['type'], '__name__'):
                    topt = '\n    type={},'.format(opt['type'].__name__)
                elif opt['type'] == 'integers':
                    topt = '\n    nargs="*",\n    type=int,'
                elif opt['type'] == 'numbers':
                    topt = '\n    nargs="*",\n    type=float,'
                elif opt['type'] == 'integer':
                    topt = '\n    type=int,'
                elif opt['type'] == 'number':
                    topt = '\n    type=float,'
                elif opt['type'] == 'string':
                    pass
                elif opt['type'] == 'strings':
                    topt = '\n    nargs="*",'
                else:
                    topt = '\n    type={},'.format(opt['type'])
            else:
                topt = '\n    nargs="*",'
        print('args.add_argument("--{}",{}{}\n    help="""{}""")'.format(name, default, topt, help))
    print('args = args.parse_args()')


if __name__ == '__main__':
    parser = argparse.ArgumentParser('simuOptParam2argparse.py',
                                     description='''This script converts a now deprecated simuOpt.Param
            list to argparse definitions. It is intended to be a
            semi-automatic process that assist your conversion from
            simuOpt.Param to argparser.
        ''')
    parser.add_argument('filename', help='''file that contains the 
        simuOpt option list.''')
    parser.add_argument('startline', type=int,
                        help='''starting line of the option list''')
    parser.add_argument('endline', type=int,
                        help='''last line of the option list''')
    parser.add_argument('name',
                        help='''Name of the option list.''')
    args = parser.parse_args()
    #
    convert_options(args.filename, args.startline, args.endline, args.name)