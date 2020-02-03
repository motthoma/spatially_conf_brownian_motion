#file that creates header which collects header for 
#confinement and intra-particle interaction which are
#used in simulation.

import sys
import os

def write_conf_int_header(conf_string, int_string):
    #function that creates header file that contains
    #includes header for confinement and inter-particle
    #interactions

    out_h = open("comp_gen_header.h", 'w')

    out_h.write('\n')

    header_interact = 'int_{}.h'.format(int_string)
    out_h.write('#include "{}"\n'.format(header_interact))

    header_conf = 'conf_{}.h'.format(conf_string)
    out_h.write('#include "{}"\n'.format(header_conf))

    out_h.close()

def call_make_file(conf_string, int_string):
    os.system("make CONF={} INT={}".format(conf_string, int_string))

if __name__ == "__main__":

    args = list(sys.argv)

    conf_string = args[1]
    int_string = args[2]

    write_conf_int_header(conf_string, int_string)
    call_make_file(conf_string, int_string)
