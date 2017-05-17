#ENO coefficients for uniform mesh
function unif_crj(k::Int)
  if k == 1
    return([1;1])
  end
  crj = zeros(k+1,k)
  for i = 1:(k+1)
    r = i-2
    for j = 1:k
      psum = 0.0
      for m = j:k
        numer = sum([*([r-q+1 for q in 0:k if q!=m && q!=l]...) for l in 0:k if l != m])
        denom = *([m-l for l in 0:k if l != m]...)
        psum = psum + numer/denom
      end
      crj[i,j] = psum
    end
  end
  crj
end

#Eno reconstruction for uniform mesh
function ENO_urec(dx,vloc::Vector,k::Int, crj::Matrix)
  vdiffs = Vector{typeof(vloc)}(0)
  vl = zero(eltype(vloc))
  vr = zero(eltype(vloc))
  N = size(vloc,1)
  if (N != 2*k-1)
    throw("dimension of vloc is not consistent with order $k ENO")
  end
  # Calculation of divided differences
  push!(vdiffs,copy(vloc))
  for i in 2:k
    push!(vdiffs,(vdiffs[i-1][2:end]-vdiffs[i-1][1:end-1])/dx)
  end
  #Calculation of Stencil
  r = 0
  if k < 2;return(vloc[1],vloc[1]);end
  for m = 2:k
    if (abs(vdiffs[m][k-r-1])<abs(vdiffs[m][k-r]))
      r = r + 1
    end
  end
  for j = 0:(k-1)
    vl = vl + vloc[k-r+j]*crj[r+1,j+1]
    vr = vr + vloc[k-r+j]*crj[r+2,j+1]
  end
  return(vl,vr)
end

#Eno reconstruction for unestructured mesh
function ENO_rec(xloc::Vector,vloc::Vector,k::Int, crj::Matrix)
  vdiffs = Vector{typeof(vloc)}(0)
  vl = zero(eltype(vloc))
  vr = zero(eltype(vloc))
  N = size(vloc,1)
  if (N != 2*k-1)
    throw("dimension of vloc is not consistent with order $k ENO")
  end
  # Calculation of divided differences
  push!(vdiffs,copy(vloc))
  for i in 2:k
    push!(vdiffs,(vdiffs[i-1][2:end]-vdiffs[i-1][1:end-1])./(xloc[i:N]-xloc[1:N-i+1]))
  end
  #Calculation of Stencil
  r = 0
  if k < 2;return(vloc[1],vloc[1]);end
  for m = 2:k
    if (abs(vdiffs[m][k-r-1])<abs(vdiffs[m][k-r]))
      r = r + 1
    end
  end
  for j = 0:(k-1)
    vl = vl + vloc[k-r+j]*crj[r+1,j+1]
    vr = vr + vloc[k-r+j]*crj[r+2,j+1]
  end
  return(vl,vr)
end
