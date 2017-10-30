class Jacobi
  require 'nmatrix'
  require 'pp'
  N, Alpha, Eps = 20, 0.04, 0.000001
  mat = NMatrix.zeros(N)
  N.times {|i| N.times {|j| N.times {|k| mat[i, j] += (Alpha * (k + 1)) ** (i + j)}}}
  mat[0, 0] = 21
  pp mat
  begin
    1.upto(N-1) do |k| k.times do |l|
        if mat[k,k] - mat[l, l] == 0
          c = s = 0.5 ** 0.5
        else
          mu = 2.0 * mat[k, l] / (mat[k, k] - mat[l, l])
          c = (0.5 * (1 + 1.0 / (1 + mu ** 2) ** 0.5)) ** 0.5
          s = mu.abs/mu * (1 - c ** 2) ** 0.5
        end
        tmp, tmp[k, k], tmp[k, l], tmp[l, k], tmp[l, l] = NMatrix.eye(N), c, -s, s, c
        mat = tmp.transpose.dot(mat).dot(tmp)
      end
    end
    ma = 0
    N.times {|i| N.times {|j| ma = mat[i, j] if mat[i, j].abs >= ma && i != j}}
  end until ma < Eps
  pp mat
end
