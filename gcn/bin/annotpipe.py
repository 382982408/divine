"""
.. module:: annotpipe
    :platform: Unix, Windows, MacOSX
    :synopsis: A wraper to call VARANT

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

This modules is a wrapper to call VARANT. The inputs are -
1. Unannotated VCF file path
2. Path to create annotated vcf file (Option)
"""
import os
import argparse
import time
import datetime
from gcn.config import lib_config
from gcn.etc.fileconfig import getlogger, FILECONFIG
from gcn.lib.varann.vartype.varant import annotator

class AnnotPline:

    def __init__(self, invcf, capture_kit_name, probe_ext_bp, outvcf=None, log_dir=None):
        self.invcf = invcf
        ts = datetime.datetime.fromtimestamp(
                            time.time()).strftime('%Y%m%d_%H%M')
        
        self.capture_kit_name =  capture_kit_name
        self.probe_ext_bp = probe_ext_bp
        

        if outvcf:
            self.outvcf = outvcf
        else:
            self.outvcf = os.path.splitext(invcf)[0] + '.varant_%s.vcf' % ts

        if log_dir:
            FILECONFIG['LOGFILE'] = os.path.join(log_dir,'varant_%s.log' % ts)
        elif outvcf:
            FILECONFIG['LOGFILE'] = os.path.splitext(self.outvcf)[0] + \
                                    '.varant_%s.log' % ts
        else:
            FILECONFIG['LOGFILE'] = os.path.splitext(invcf)[0] +\
                                    '.varant.%s.log' % ts

        self.logger = getlogger()
        status = self._check_env()  # Check Environment variable setting
        if status == 1:
            d = 'There seems to be problem with setting the '\
                'environment variables..'
            self.logger.error(d)
            print d
        else:
            self.logger.info('Environment variable check successful..')

    def _check_env(self):
        flag = 0
        envs = ['GCN', 'GCN_DATA_DIR', 'GCN_DB_DIR']
        for env in envs:
            loc = lib_config.gcn_path(env)
            if loc:
                if not os.path.exists(loc):
                    self.logger.error('Path set to %s environment does not '
                                      'exist..' % env)
                    flag = 1
            else:
                self.logger.error('%s environment variable is not set..' % env)
                flag = 1
        return flag

    def annotate_varant(self,hgmd_on=False):
        'Annotate the VCF using VARANT'
        self.logger.info('Input file = %s, Output file = %s' %
                         (self.invcf, self.outvcf))
        annotator.main(self.invcf, self.outvcf, self.capture_kit_name, self.probe_ext_bp, hgmd_on, self.logger)
        d = 'Annotation complete [%s;%s].'%(self.invcf,self.outvcf)
        print d
        self.logger.info(d)

def main():
    """Main script to annotate VCF file"""

    desc = 'VCF Annotator'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--inputVCF', dest='invcf', type=str,
                        help='Input VCF file path')
    parser.add_argument('-c', '--captureKitName', dest='capture_kit_name', type=str, default = 'SeqCapEZ_Exome',
                        help='Input Capture kit name')
    parser.add_argument('-e', '--probeFlankingB', dest='probe_flanking_bp', type=int, default = 50,
                        help='Input Capture kit probe extension bp')
    parser.add_argument('--hgmd', action='store_const', dest='hgmd', required=False, default=False, const=True, help='enable HGMD (requires a license)[False]')
    parser.add_argument('-o', '--outVCF', dest='outvcf', type=str,
                        default=None, help='Output VCF file path')
    parser.add_argument('-l', '--logDir', dest='logdir', type=str,
                        default=None, help='Log directory file path')
    args = parser.parse_args()
    ap = AnnotPline(args.invcf, args.capture_kit_name, args.probe_flanking_bp, args.outvcf, args.logdir)
    
    ap.annotate_varant(args.hgmd)


if __name__ == "__main__":
    main()
