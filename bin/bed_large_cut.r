#! /usr/bin/env Rscript

# needlestack: a multi-sample somatic variant caller
# Copyright (C) 2015 Matthieu Foll and Tiffany Delhomme
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
cut_cum_size=floor(seq(1,bed_cum_size[length(bed_cum_size)],l=nb_pieces+1))[-c(1)]
lines_cut_bed=unique(c(0,unlist(lapply(cut_cum_size,function(x){min(which(bed_cum_size>=x))}))))
for (k in 1:(length(lines_cut_bed)-1)) {
  cur_region=(lines_cut_bed[k]+1):lines_cut_bed[k+1]
  region_name=paste(bed[cur_region[1],1],"_",bed[cur_region[1],2],"-",bed[cur_region[length(cur_region)],1],"_",bed[cur_region[length(cur_region)],3],sep="")
  if (is.na(commandArgs(TRUE)[2])) {
    cat(paste(bed[cur_region,1],":",bed[cur_region,2],"-",bed[cur_region,3],sep=""),sep="\n",file=paste(region_name,"_regions",sep=""))
  } else {
  	cat(paste(bed[cur_region,1],"\t",bed[cur_region,2],"\t",bed[cur_region,3],sep=""),sep="\n",file=paste(region_name,".bed",sep=""))
  }
}
