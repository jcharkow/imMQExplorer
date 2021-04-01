# The following functions are not specific to an analyses but used in several scripts to make code simpler to read

#convert string argument to bool (using code from stack overflow) https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
import argparse
def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


