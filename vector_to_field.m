function [field] = vector_to_field(vector,i,j,k)

NZ = max(k);
NY = max(j);
NX = max(i);
NVEC = length(vector);
field = nan(NZ,NY,NX);
for nv = 1:NVEC
 field(k(nv),j(nv),i(nv)) = vector(nv);
end
