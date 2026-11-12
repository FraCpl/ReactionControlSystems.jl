"""
    My = rcsMixMatrix(posRCS, dirRCS, forceRCS, posRef)

Compute the RCS mixing matrix ```My``` with respect to ```posRef```, so that:

[Force; Torque] = ```My```*y

with 0 ≤ y ≤ 1.
"""
function rcsMixMatrix(posRCS, dirRCS, forceRCS=ones(length(posRCS)), posRef=zeros(3))
    n = length(posRCS)
    dir = normalize.(dirRCS)
    pos = posRCS .- Ref(posRef)
    My = zeros(6, n)
    if length(forceRCS) == 0
        frc = ones(length(posRCS))
    else
        frc = forceRCS
    end
    for i in 1:n
        My[:, i] = frc[i]*[dir[i]; pos[i] × dir[i]]
    end
    return My
end
