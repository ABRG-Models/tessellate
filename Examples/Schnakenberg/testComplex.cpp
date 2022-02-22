using namespace std::complex_literals;
    Vector<std::complex<double>, 4> cplx;
    cplx.set_from (std::pow(1i, 2));
    std::cout << "Complex*2: " << cplx*2.0 << std::endl;
