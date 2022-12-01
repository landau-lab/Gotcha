#' Merge multi-sample Gotcha fastqs
#' @param path Path to Gotcha fastqs
#' @param out.dir Name of new directory to save merged fastqs
#' @return New directory with merged fastqs
#'

MergeFastqs = function(path='gotcha_fastqs/', out.dir='merged'){

	setwd(path)

	bash_command = 'for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)

		do echo "Merging R1"

	cat "$i"_L00*_R1_001.fastq.gz > merged_"$i"_L001_R1_001.fastq.gz

		echo "Merging R2"

	cat "$i"_L00*_R2_001.fastq.gz > merged_"$i"_L001_R2_001.fastq.gz

		echo "Merging R3"

	cat "$i"_L00*_R3_001.fastq.gz > merged_"$i"_L001_R3_001.fastq.gz

	done;'

	message("------- START MERGE FASTQS -------")

	system(bash_command)

	message("------- DONE MERGING -------")

	system('mkdir merged')

	message("------- MOVING MERGED FASTQS TO NEW FOLDER -------")

	system(paste0('mv merged_*', out.dir))

	}

