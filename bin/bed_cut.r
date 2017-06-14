#! /usr/bin/env Rscript

# needlestack: a multi-sample somatic variant caller
# Copyright (C) 2015 IARC/WHO
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

nb_pieces=as.numeric(commandArgs(TRUE)[1])
bed=read.table(pipe('cat /dev/stdin'),stringsAsFactors = F)
bed_cum_size=cumsum(bed[,3]-bed[,2]+1)
if(nb_pieces>bed_cum_size[length(bed_cum_size)]) nb_pieces=bed_cum_size[length(bed_cum_size)]
if(nb_pieces==1) { cat(paste(bed[,1],":",bed[,2],"-",bed[,3], sep=""), sep="\n", file=paste(bed[1,1],"_",bed[1,2],"-",bed[dim(bed)[1],1],"_",bed[dim(bed)[1],3],"_regions",sep="")) } else {
  cut_cum_size=floor(seq(1,bed_cum_size[length(bed_cum_size)],l=nb_pieces+1))
  for(k in 1:(nb_pieces-1) ) {
    cur_size = cut_cum_size[k+1]
    cut_id = which(bed_cum_size > cur_size)[1]
    if (length(bed_cum_size)==1 | cut_id == 1) {new_cut = bed[cut_id,2] + cur_size} else { new_cut = bed[cut_id,2] + (cur_size - bed_cum_size[cut_id-1]) }
    if(k==1) {cut_data=data.frame(c(cut_id,new_cut))} else {cut_data = cbind(cut_data,c(cut_id,new_cut))}
  }
  for(i in 1:nb_pieces ){
    if (i==1){
      chrs=bed[1:cut_data[1,i],1]
      starts=bed[1:cut_data[1,i],2]
      ends=bed[1:cut_data[1,i],3]
      ends[length(ends)]=cut_data[2,i]-1
      if(cut_data[2,i]-1<starts[length(starts)]) { chrs=chrs[1:(length(chrs)-1)];starts=starts[1:(length(starts)-1)];ends=ends[1:(length(ends)-1)] }
      cat(paste(chrs,":",starts,"-",ends, sep=""), sep="\n", file=paste(chrs[1],"_",starts[1],"-",chrs[length(chrs)],"_",ends[length(ends)],"_regions",sep=""))
    } else if (i==nb_pieces) {
      chrs=bed[cut_data[1,i-1]:dim(bed)[1],1]
      starts=bed[cut_data[1,i-1]:dim(bed)[1],2]
      starts[1]=cut_data[2,i-1]
      ends=bed[cut_data[1,i-1]:dim(bed)[1],3]
      cat(paste(chrs,":",starts,"-",ends, sep=""), sep="\n", file=paste(chrs[1],"_",starts[1],"-",chrs[length(chrs)],"_",ends[length(ends)],"_regions",sep=""))
    } else {
      chrs=bed[cut_data[1,i-1]:cut_data[1,i],1]
      starts=bed[cut_data[1,i-1]:cut_data[1,i],2]
      ends=bed[cut_data[1,i-1]:cut_data[1,i],3]
      starts[1]=cut_data[2,i-1]
      ends[length(ends)]=cut_data[2,i]-1
      if(cut_data[2,i]-1<starts[length(starts)]){ chrs=chrs[1:(length(chrs)-1)];starts=starts[1:(length(starts)-1)];ends=ends[1:(length(ends)-1)] }
      cat(paste(chrs,":",starts,"-",ends, sep=""), sep="\n", file=paste(chrs[1],"_",starts[1],"-",chrs[length(chrs)],"_",ends[length(ends)],"_regions",sep=""))
    }
  }
}
