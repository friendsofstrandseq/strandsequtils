'''
This module converts genome STRIP CNV VCF format to a bed file with some extra columns (for SVs)
'''

import sys

def output_bed_format_CNV(vcf_file, bed_file, sample):
	
	'''
	This function reads the genome STRIP CNV vcf file, and outputs the SV calls in a 6-column format (chrom, start, end, SV type, GT, CN) for the required sample
	
	Args:
		vcf_file: the input vcf file containing CNVs
		bed_file: the output file in bed format
		sample: the sample name to be processed in the VCF file
	'''
	
	with open(vcf_file) as f:
		with open(bed_file, 'w') as out:
			out.write('chrom\tstart\tend\tSV\tGT\tCN\n')
			for line in f:
				# skip the header lines
				if line[0]=='#':
					columns = line
					continue
					
				sp_columns = columns.split()
							
				if not sample in sp_columns:
					print('Error:', sample, 'is not present in the header columns')
				else:
					sample_col_index = sp_columns.index(sample)
					
				# processing non header lines
				sp = line.split()
				sv_type = ''
				end = '-1'
				
				chrom = sp[0]
				start = sp[1] 
				info = sp[7]
				format = sp[8]
				
				info_sp = info.split(';')
				
				# finding END pos and SV TYPE
				for s in info_sp:
					if s.startswith('SVTYPE='):
						sv_type = s[7:]
					elif s.startswith('END='):
						end = s[4:]
						
				# finding GT and Cn indiced in the colon separated fields
				format_sp = format.split(':')
				
				gt = format_sp.index('GT') if 'GT' in format_sp else None
				cn = format_sp.index('CN') if 'CN' in format_sp else None
				
							
				# processing the sample column and getting GT and CN values
				sample_col = sp[sample_col_index]
				sample_col_sp = sample_col.split(':')
				
				sample_gt = sample_col_sp[gt] if gt is not None else '.'
				sample_cn = sample_col_sp[cn] if cn is not None else '.'

				if sv_type=='CNV' and sample_cn=='2':
					# there is no CNV in this sample
					continue
				
				out.write(chrom + '\t' + start + '\t' + end + '\t' + sv_type + '\t' + sample_gt + '\t'  + sample_cn + '\n')
				
					
					

def main(args):
	output_bed_format_CNV(args[0], args[1], args[2])


if __name__=='__main__':
	main(sys.argv[1:])
