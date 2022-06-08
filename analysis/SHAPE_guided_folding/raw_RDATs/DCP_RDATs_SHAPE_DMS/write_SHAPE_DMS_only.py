from glob import glob
import os, sys

files=glob('../RDAT/*.rdat')

for fil in files:
	with open(os.path.basename(fil),'w') as outfile:
		with open(fil,'r') as infile:
			for lin in infile.readlines():
				if lin.startswith('ANNOTATION_DATA:') or lin.startswith('REACTIVITY:') or lin.startswith('REACTIVITY_ERROR'):
					if ':9' in lin or ':12' in lin:
						#print('this1', lin)
						outfile.write(lin.replace(':9',':1').replace(':12',':2'))

				else:
					#print('this2', lin)
					outfile.write(lin)
