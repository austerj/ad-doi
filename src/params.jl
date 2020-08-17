# parameter set H
s₀ = 100
ν₀ = 0.16
r = 0.04
κ = 0.6
θ = 0.04
ξ = 0.2
ρ = -0.15
heston = Heston(s₀, ν₀, r, κ, θ, ξ, ρ)

# multivariate parameter set H
m_s₀ = [100; 90]
m_ν₀ = [0.02; 0.16]
m_κ = [0.6; 0.6]
m_θ = [0.04; 0.04]
m_ξ = [0.2; 0.2]
m_ρ = [1    -0.15 -0.4  0.25;
      -0.15  1     0.3 -0.2;
      -0.4   0.3   1    0.15;
       0.25  -0.2  0.15 1]
m_heston = MultivariateHeston(m_s₀, m_ν₀, r, m_κ, m_θ, m_ξ, m_ρ)
