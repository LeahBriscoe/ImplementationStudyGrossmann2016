import csv
import fnmatch, os, re
import numpy as np



in_paths = []
short_in_paths = []
ref_path = 'GrossmannsModulescsv'


module_ref = {}
my_modules = {}


for file in os.listdir('../ISA_Modules_Results/'):
	if fnmatch.fnmatch(file, 'Module*.txt'):
		short_in_paths.append(str.replace(file,"/",''))
		in_paths.append(os.path.join('../ISA_Modules_Results/', file))

out_paths = [str('OverlapOutput'+ in_x) for in_x in short_in_paths]

with open(ref_path, 'r') as csvfile:
	my_reader = csv.reader(csvfile)
	for row in my_reader:
		row = filter(None,row)
		splitpoint = row.index("ImagingFeature")
		module_ref[row[0]] = {}
		module_ref[row[0]]['Pathway'] =  row[1:splitpoint]
		module_ref[row[0]]['ImagingFeature'] = row[splitpoint+1:len(row)]

for p in range(0,len(in_paths)):
	with open(in_paths[p],'r') as my_module_file:
		my_module_lines= my_module_file.read().splitlines()
	my_sep = my_module_lines.index('')
	my_modules[p] = {}
	my_modules[p]['Pathway'] = my_module_lines[0:my_sep]
	my_modules[p]['ImagingFeature'] = my_module_lines[(my_sep+2):len(my_module_lines)]

	clean1_pathways = [str.replace(pa,"REACTOME_",'') for pa in my_modules[p]['Pathway']]
	clean_pathways = [str.replace(pa,"_",' ') for pa in clean1_pathways]


	print in_paths[p]
	# problem is alphabetical sorting doesn't match!

	sorted_keys = ["M1_module1","M2_module2","M3_module3","M4_module4","M5_module5","M6_module6","M7_module7",
		"M8_module8","M9_module9","M10_module10","M11_module11","M12_module12","M13_module13"]
	pathway_similarity = []
	pathway_sizes = []
	pathway_commonality = []

	imagingfeature_similarity = []
	imagingfeature_sizes = []
	imagingfeature_commonality= []
	for Mx in sorted_keys:
		
		print Mx
		pathway_similarity.append(float(len(set(clean_pathways).intersection(module_ref[Mx]['Pathway'])))/len(my_modules[p]['Pathway']))
		pathway_sizes.append(tuple((len(set(clean_pathways).intersection(module_ref[Mx]['Pathway'])),len(module_ref[Mx]['Pathway']))))
		pathway_commonality.append(set(clean_pathways).intersection(module_ref[Mx]['Pathway']))

		imagingfeature_similarity.append(float(len(set(my_modules[p]['ImagingFeature']).intersection(module_ref[Mx]['ImagingFeature'])))/len(my_modules[p]['ImagingFeature']))
		imagingfeature_sizes.append(tuple((len(set(my_modules[p]['ImagingFeature']).intersection(module_ref[Mx]['ImagingFeature'])),len(module_ref[Mx]['ImagingFeature']))))

		imagingfeature_commonality.append(set(my_modules[p]['ImagingFeature']).intersection(module_ref[Mx]['ImagingFeature']))

	print max(pathway_similarity)
	print pathway_similarity
	print max(imagingfeature_similarity)
	print imagingfeature_similarity
	max_path_index = pathway_similarity.index(max(pathway_similarity))
	max_feature_index = imagingfeature_similarity.index(max(imagingfeature_similarity))

	with open(out_paths[p], 'w') as fout:
		fout.write(str("Max Pathway Similarity " + str(int(max_path_index+1)) + " at " + str(100* max(pathway_similarity)) + "%" + 
			"\n" + "Max Imaging Feature Similarity " + str(int(max_feature_index+1)) + " at " + str(100*max(imagingfeature_similarity)) + "%"))
		fout.write("\nNo. Pathways in common and Orig Module Size" + str(pathway_sizes[max_path_index]) +"\n" + "No.  Imaging Features in common and Orig module size" + str(imagingfeature_sizes[max_feature_index]))
		fout.write("\nPathways in common")
		for i in pathway_commonality[max_path_index]:
			fout.write("\n")
			fout.write(i)
		fout.write("\n \nImaging Features in common")
		for j in imagingfeature_commonality[max_feature_index]:
			fout.write("\n")
			fout.write(j)
		
		fout.close()




		
