import numpy as np
import argparse
import time

def write_nve(outf,pot_dir,pot_name,index,datainp,temp):
    with open(outf, 'w') as of :
        of.writelines("log log{}.log".format(index))
        of.writelines("units metal")
        of.writelines("atom_style atomic")
        of.writelines("read_data  {}".format(datainp))
        of.writelines("pair_style deepmd {}/{}".format(pot_dir,pot_name))
        of.writelines("pair_coeff")

        if dis == True:
            randomnum = time.time_ns()
            of.write("displace_atoms all random 0.02 0.02 0.02 {} units box".format(randomnum))
            of.write("velocity all create {} {} mom yes rot yes dist gaussian".format(temp,randomnum))

        of.writelines("timestep {}e-03".format(dt))
        of.writelines("thermo {}".format(freq_t))
        of.writelines("thermo_style custom step time ke pe etotal temp press density lx ly lz")
        of.writelines("fix finvt1 all nvt temp {} {} $(100.0*dt)".format(temp,temp))

        if dump == True:
            of.writelines("dump NVTdump all custom {} dump{}.lammpstrj id type x y z ix iy iz vx vy vz fx fy fz mass".format(freq_d*10,index))
            of.writelines("dump binDump4 all custom {} dump{}.bin id type x y z vx vy vz".format(freq_d,index))
            of.writelines("run {}".format(numstep))

        of.writelines("write_restart final{}.restart".format(index))
        of.writelines("write_data final{}.data".format(index))
    return
def write_nvt():
    return
if(__name__ == "main"):
    parser = argparse.ArgumentParser(description='Generate NVE input files')
    parser.add_argument('-d', '--directory',
                        type = ,
                        required = ,
                        help = ,
                        default = )