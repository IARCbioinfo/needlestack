region="17:7571720-7579962"
samples=readLines(file("CEU_samples.txt"))
for (cur_sample in samples) {
  print(cur_sample)
  cur_filename=readLines(pipe(paste("curl -l ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/",cur_sample,"/exome_alignment/ --user anonymous:follm@iarc.fr | grep .*\\\\.mapped.*bam$",sep="")))
  if (length(cur_filename)>0) {
    system(paste("samtools view -b -h -o ",cur_sample,".bam ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/",cur_sample,"/exome_alignment/",cur_filename," 17:7571720-7579962",sep=""))
    system(paste("samtools index ",cur_sample,".bam",sep=""))
    system("rm *mapped*.bai")
  } else {
    print("...missing")
  }
}