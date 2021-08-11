function [vector] = field_to_vector(field,i,j,k)

  Nfield = length(i);
  vector = zeros(Nfield,1);
  for ni = 1:Nfield
    vector(ni) = field(k(ni),j(ni),i(ni));
  end

