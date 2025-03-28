export gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, gamma9


function gamma1(d, jmm1)
    j, m, m1 = jmm1
    N = d.N * 1.0

    spontaneous = d.collective_emission / 2 * (2 * j * (j + 1) - m * (m - 1) - m1 * (m1 - 1))
    losses = (d.emission / 2) * (N + m + m1)
    pump = d.pumping / 2 * (N - m - m1)
    collective_pump = d.collective_pumping / 2 * (2 * j * (j + 1) - m * (m + 1) - m1 * (m1 + 1))
    collective_dephase = d.collective_dephasing / 2 * (m - m1)^2

    if j <= 0
        dephase = d.dephasing * N / 4
    else
        dephase = d.dephasing / 2 * (N / 2 - m * m1 * (N / 2 + 1) / j / (j + 1))
    end

    return -(spontaneous + losses + pump + dephase + collective_pump + collective_dephase)

end

function gamma2(d::Dicke, jmm1::Tuple{Number,Number,Number})
    j, m, m1 = jmm1
    N = d.N * 1.0

    if d.collective_emission == 0
        spontaneous = 0.0
    else
        spontaneous = d.collective_emission * sqrt((j + m) * (j - m + 1) * (j + m1) * (j - m1 + 1))
    end

    if ((d.emission == 0) || (j <= 0))
        losses = 0.0
    else
        losses = d.emission / 2 * sqrt((j + m) * (j - m + 1) * (j + m1) * (j - m1 + 1)) * (N / 2 + 1) / (j * (j + 1))
    end

    return spontaneous + losses

end

function gamma3(d::Dicke, jmm1::Tuple{Number,Number,Number})
    j, m, m1 = jmm1
    N = d.N * 1.0

    if ((d.emission == 0) || (j <= 0))
        return 0.0
    else
        return d.emission / 2 * sqrt((j + m) * (j + m - 1) * (j + m1) * (j + m1 - 1)) * (N / 2 + j + 1) / (j * (2 * j + 1))
    end
end

function gamma4(d::Dicke, jmm1::Tuple{Number,Number,Number})
    N = d.N * 1.0
    j, m, m1 = jmm1
    if ((d.emission == 0) || ((j + 1) <= 0))
        return 0.0
    else
        return d.emission / 2 * sqrt((j - m + 1) * (j - m + 2) * (j - m1 + 1) * (j - m1 + 2)) * (N / 2 - j) / ((j + 1) * (2 * j + 1))
    end
end


function gamma5(d::Dicke, jmm1::Tuple{Number,Number,Number})
    N = d.N * 1.0
    j, m, m1 = jmm1
    if ((d.dephasing == 0) || (j <= 0))
        return 0.0
    else
        return d.dephasing / 2 * sqrt((j^2 - m^2) * (j^2 - m1^2)) * (N / 2 + j + 1) / (j * (2 * j + 1))
    end
end


function gamma6(d::Dicke, jmm1::Tuple{Number,Number,Number})
    j, m, m1 = jmm1
    N = d.N * 1.0
    if d.dephasing == 0
        return 0.0
    else
        return d.dephasing / 2 * sqrt(((j + 1)^2 - m^2) * ((j + 1)^2 - m1^2)) * (N / 2 - j) / ((j + 1) * (2 * j + 1))
    end
end

function gamma7(d::Dicke, jmm1::Tuple{Number,Number,Number})
    j, m, m1 = jmm1
    N = 1.0 * d.N
    if ((d.pumping == 0) || (j <= 0))
        return 0.0
    else
        return d.pumping / 2 * sqrt((j - m - 1) * (j - m) * (j - m1 - 1) * (j - m1)) * (N / 2 + j + 1) / (j * (2 * j + 1))
    end
end


function gamma8(d::Dicke, jmm1::Tuple{Number,Number,Number})
    j, m, m1 = jmm1
    N = 1.0 * d.N

    if ((d.pumping == 0) || (j <= 0))
        pump = 0.0
    else
        pump = d.pumping / 2 * sqrt((j + m + 1) * (j - m) * (j + m1 + 1) * (j - m1)) * (N / 2 + 1) / (j * (j + 1))
    end

    if (d.collective_pumping == 0)
        collective_pump = 0.0
    else
        collective_pump = d.collective_pumping * sqrt((j - m) * (j + m + 1) * (j + m1 + 1) * (j - m1))
    end
    return pump + collective_pump
end

function gamma9(d::Dicke, jmm1::Tuple{Number,Number,Number})
    j, m, m1 = jmm1
    N = d.N * 1.0
    if (d.pumping == 0)
        return 0.0
    else
        return d.pumping / 2 * sqrt((j + m + 1) * (j + m + 2) * (j + m1 + 1) * (j + m1 + 2)) * (N / 2 - j) / ((j + 1) * (2 * j + 1))
    end
end