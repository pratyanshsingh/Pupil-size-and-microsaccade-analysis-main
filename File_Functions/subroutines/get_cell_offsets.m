function [cell_data_offsets]=get_cell_offsets(fhdr, chdr)
%function [cell_data_offsets]=get_cell_offsets(fhdr, chdr)
%
%Returns the index to the beginning of each cell's data
%

%Check HdrVer to see if file is a ses file or ave file
if fhdr(2) == 1
  ses_hdr_offsets_v;
else
  ave_hdr_offsets_v;
end

cell_data_offsets = zeros(1,fhdr(NCells));
cum_cell_sizes = 0;

for c = 1:fhdr(NCells)
	cell_data_offsets(c) = fhdr(LHeader) + cum_cell_sizes;
	cum_cell_sizes = cum_cell_sizes + chdr(c,NObs)*chdr(c,NPoints)*fhdr(NChan)*2;
end
