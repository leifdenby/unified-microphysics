import glob
import re
import argparse
import datetime

WRAPPER_MODULE_NAME = 'microphysics_initialisation'
CONFIGURATION_NAMELIST_FILENAME = 'input/microphysics.nml'

parser = argparse.ArgumentParser('Fortran wrapper compiler for unified microphysics')
parser.add_argument('--debug', default=False, action='store_true')
args = parser.parse_args()

debug = args.debug

files = glob.glob('mphys_*.F*')
mphys_names = [re.match(r'mphys_(?P<mphys_name>.*).F[\d2]', file).group('mphys_name') for file in files]

if debug:
    print ", ".join(mphys_names)


print "Found %d microphysics implementations. Building wrapper..." % len(files)


MODULE_TEMPLATE = """! Generated by unified microphysics wrapper generator on {datetime}

""" + open('templates/%s.tpl.F90' % WRAPPER_MODULE_NAME).read().replace('!{', '{')

MODULE_SUBROUTINE_TEMPLATE = """
    subroutine mphys_wrapper_{mphys_name}
        use mphys_{mphys_name}, only: init, calc_dq !, n_moments
        !n_tracers = size(n_moments)
        !if (maxval(n_moments) > n_moments__max) then
        !    n_moments__max = maxval(n_moments)
        !endif
        !call allocate_local()
        call init()
        q_flux_function => calc_dq
    end subroutine mphys_wrapper_{mphys_name}
"""

SWITCH_STATEMENT_TEMPLATE = """
    else if (configuration .eq. "{mphys_name}") then
        call mphys_wrapper_{mphys_name}()
    """


subroutine_blocks = '\n'.join([MODULE_SUBROUTINE_TEMPLATE.format(mphys_name=m) for m in mphys_names])

template_kwargs = {
    'datetime': datetime.datetime.now().isoformat(),
    'wrapper_name': WRAPPER_MODULE_NAME,
    'microphysics_configuration_namelist_filename': CONFIGURATION_NAMELIST_FILENAME,
    'microphysics_module_blocks': subroutine_blocks,
    'mphys_names_joined': ', '.join(mphys_names),
    'mphys_methods_switch_statements': '\n'.join([SWITCH_STATEMENT_TEMPLATE.format(mphys_name=m) for m in mphys_names]),
}

output_filename = '%s.F90' % WRAPPER_MODULE_NAME
with open(output_filename, 'w') as f:
    f.write(MODULE_TEMPLATE.format(**template_kwargs))

if debug:
    print "Wrapper written to %s" % output_filename



