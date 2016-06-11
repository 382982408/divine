#!/usr/bin/env python
'''
author: changjin.hong@gmail.com
'''

import os, sys, urllib, requests, argparse, shutil, imp, re
import subprocess as sp

VERSION = '0.1.1'
author_email = 'changjin.hong@gmail.com'

def cmd_exists(cmd):
	return sp.call("type " + cmd, shell=True,\
		stdout=sp.PIPE, stderr=sp.PIPE) == 0

def has_module(module):
	try:
		imp.find_module(module)
		print '[%s] is in path or installed'%module
		found = True
	except ImportError:
		print '[%s] is not in path or installed'%module
		found = False
	return found
	
def syscmd(cmd):
	subproc = sp.Popen(cmd, shell=True)
	output,error = subproc.communicate()
	exit_code = subproc.wait()
	if exit_code>0:
		raise RuntimeError('command[%s] failed;%s'%(cmd,error))
	
def check_py_ver():
	verStr = '%d.%d'%(sys.version_info[0],sys.version_info[1])
	return verStr

def check_network_on():
	try :
		stri = "https://www.google.com"
		data = urllib.urlopen(stri)
		print "Connected"
		return True
	except e:
		print "not connected" ,e 
		return False

def get_pylib_prefix():
	return 'python%s'%check_py_ver()

class Setup():
	def __init__(self,url_fn=None):
		#to get the current directory
		self.install_dir = os.path.dirname(os.path.abspath(__file__))
		self.py_libs_dir = os.path.join(self.install_dir,'python_libs','pkg')
		self.py_libs_prefix = ['fastsemsim-0.9.4','hgvs','hpo_similarity']
		self.py_libs = ['fastsemsim','pyhgvs','hpo_similarity']
		
		self.py_modules = ['ConfigParser','backports-abc',\
			'backports.ssl-match-hostname','certifi','decorator',\
				'matplotlib','networkx','nose','numpy','pandas','pygr',\
					'pyparsing','pysam','python-dateutil','pytz',\
						'scikit-learn','singledispatch','six','tornado','xlwt']

		self.pylib_prefix = 'python_libs'
		self.pylib = os.path.join(self.install_dir,self.pylib_prefix)
		self.pylib_build_prefix = 'site-packages'
		self.pylib_build = os.path.join(self.pylib,'lib',get_pylib_prefix(),self.pylib_build_prefix)
		
		if not os.path.exists(self.pylib_build):
			os.makedirs(self.pylib_build)

		self.gcn = 'gcn'
		self.gcndata = 'gcndata'
		self.gcndb = 'gcndb'
		self.gcnlog = 'gcn/logs/gcn.log'
		
		if url_fn is None:
			self.url_fn = os.path.join(self.install_dir,'managment_script','support_pkg_urls.txt')
		else:
			self.url_fn = url_fn
			
		if not os.path.exists(self.url_fn):
			raise RuntimeError('a config file[%s] for resource packages does not exist!'%self.url_fn)

	def install_module_by_pip(self,pkg_name):
		print 'installing %s ...'%pkg_name
		cmd = "PYTHONUSERBASE=%s pip install --user %s"%(self.pylib,pkg_name)
		syscmd(cmd)
	
	def download_data(self,keep_download_files=False):
		print "updating prebuilt annotation databases (be patient) ..."
		
		#to access a file containing url file
		print "downloading Divine resource files ..."
		fp = open(url_fn)
		for i in fp:
			if i.startswith('#'):continue
			res_name,res_url = i.strip().split('\t')
			if cmd_exists('wget'):
				cmd = 'cd %s;wget -c %s'%(self.install_dir,res_url)
			elif cmd_exists('rsync'):
				cmd = "rsync -ravzP %s %s/"%(res_url,self.install_dir)
			else:
				raise RuntimeError('either rsync or wget should be in path!')
			syscmd(cmd)
		fp.close()
		print "Done."
		
		print "extracting resource files ..."
		fns = os.listdir(self.install_dir)
		extracted = []
		for fn in fns:
			if 'tar.bz2' in fn:
				mObj = re.search(r'(\S+)\.tar\.bz2\.\S+',fn)
				if mObj:
					fbase='%s.tar.bz2'%mObj.group(1)
					arch_path = os.path.join(self.install_dir,fbase)
					cmd = "cat %s.* | tar -xvjf - -C %s"%(arch_path,self.install_dir)
				else:
					arch_path = os.path.join(self.install_dir,fn)
					cmd = "tar -xvjf %s -C %s"%(arch_path,self.install_dir)

				if arch_path not in extracted:
					extracted.append(arch_path)
					syscmd(cmd)
		print "Done."
		if not keep_download_files:
			os.system("rm -rf %s/*.tar.bz2*"%self.install_dir)
	
	def uninstall_resource(self):
		shutil.rmtree(os.path.join(self.install_dir,'gcndb'))
		shutil.rmtree(os.path.join(self.install_dir,'gcndata'))
		
	def install_python_libs(self):
		#to make sure that PYTHONPATH includes the directory to install python modules
		if 'PYTHONPATH' in os.environ:
			for module in self.py_libs_prefix:
				module_path = os.path.join(self.py_libs_dir,module)
				setup_py = os.path.join(module_path,'setup.py')
				if os.path.exists(module_path) and os.path.exists(setup_py):
					cmd = "export PYTHONPATH=%s:%s:$PYTHONPATH"%(self.pylib,self.pylib_build)
					cmd += ";cd %s;python ./setup.py install --prefix=%s"%(module_path,self.pylib)
					syscmd(cmd)
				else:
					raise IOError('check if [%s] contains all necessary files'%module_path)
				
			for module in self.py_modules:
				if not has_module(module):
					self.install_module_by_pip(module)
		else:
			raise EnvironmentError("add PYTHONPATH into your shell configuration file.")
		
		print 'python modules are installed successfully in [%s]'%self.pylib
	
	def uninstall_python_libs(self):
		if 'PYTHONPATH' in os.environ:
			if self.pylib in os.environ['PYTHONPATH']:
				for module in self.py_libs:
					if has_module(module):
						cmd = "pip uninstall -y %s"%module
						try:
							syscmd(cmd)
						except:
							print "%s may not be installed previously."%module
					else:
						print "%s may not be installed previously."%module
				path2del = os.path.join(self.pylib,'lib')
				if os.path.exists(path2del):
					shutil.rmtree(path2del)
					
				path2del = os.path.join(self.pylib,'bin')
				if os.path.exists(path2del):
					shutil.rmtree(path2del)
			else:
				raise EnvironmentError("your PYTHONPATH does not include '%s'."%self.pylib)
		else:
			raise EnvironmentError("add PYTHONPATH into your shell configuration file.")
		
		print 'python modules are uninstalled successfully in [%s]'%self.pylib
	
	def msg_config(self,setup_mode='install'):
		msg = "export DIVINE=%s\n"%self.install_dir
		msg += "export GCN=$DIVINE/%s\n"%self.gcn
		msg += "export GCN_DATA_DIR=$DIVINE/%s\n"%self.gcndata
		msg += "export GCN_DB_DIR=$DIVINE/%s\n"%self.gcndb
		msg += "export GCN_LOGFILE=$DIVINE/%s\n"%self.gcnlog
		
		hline = "-------------------------------"

		if setup_mode == 'install':
			desc1='add'
			desc2='to'
		
		elif setup_mode == 'remove':
			desc1='remove'
			desc2='from'
		
		print "\n"
		print '-First, %s the following variables %s your shell script (e.g., $HOME/.bash_profile or $HOME/.profile).'%(desc1,desc2)
		print hline
		print msg
		
		print "-Then, %s the following variables %s PYTHONPATH in your shell script (e.g., $HOME/.bash_profile or $HOME/.profile)."%(desc1,desc2)
		print hline
		msg = "export PATH=$DIVINE:$PATH\n"
		msg += "export PYTHONPATH=$DIVINE:$PYTHONPATH\n"
		msg += "export PYTHONPATH=$DIVINE/%s:$PYTHONPATH\n"%self.pylib_prefix
		msg += "export PYTHONPATH=$DIVINE/%s/%s/%s/%s:$PYTHONPATH\n"%(self.pylib_prefix,'lib',get_pylib_prefix(),self.pylib_build_prefix)
		print msg
		
def main():
	parser = argparse.ArgumentParser(description="Divine setup (v%s) [author:%s]"%(VERSION,author_email))
	parser.add_argument('--install', action='store_const', dest='install', required=False, default=False, const=True, help='install Divine')
	parser.add_argument('--update_db', action='store_const', dest='update_db', required=False, default=False, const=True, help='update_db [False]')
	parser.add_argument('-c', action='store', dest='url_fn', required=False, default=None, help='a file containing URL to download divine resource files')
	parser.add_argument('--uninstall', action='store_const', dest='uninstall', required=False, default=False, const=True, help='uninstall Divine')
	parser.add_argument('--remove_db', action='store_const', dest='remove_db', required=False, default=False, const=True, help='remove_db [False]')

	args = parser.parse_args()

	if (args.install and args.uninstall) or (not args.install and not args.uninstall):
		raise RuntimeError('choose --install or --uninstall')
	
	#to install python libraries
	cs = Setup(url_fn=args.url_fn)
	if args.install:
		if not check_network_on():
			raise RuntimeError('it requires network connection to setup Divine!')
		cs.install_python_libs()
		if args.update_db:
			cs.download_data()
		cs.msg_config('install')
	elif args.uninstall:
		cs.uninstall_python_libs()
		if args.remove_db:
			cs.delete_resource()
		cs.msg_config('remove')
	
if __name__ == '__main__':
	main()
