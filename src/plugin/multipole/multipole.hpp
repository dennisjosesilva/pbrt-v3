class DipoleSolver
{
public:
  DipoleSolver(Float sigma_a, Float sigma_s_prime, Floaat d);

  Float R(Float r);
  Float T(Float r);
private:
  Float sigma_a;
  Float sigma_s_prime;
  Float d;
};
