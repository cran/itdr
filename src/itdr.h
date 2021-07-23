void ITM_mean(const double* x, const double* gk, const double* y,
              const double* h, const int* n, const int* p,
              const int* out, double* Mmat);

void ITM_pdf(const double* x, const double* gk, const double* y,
             const double* hx, const double* hy,
             const int* n, const int* p,
             const int* out, double* Mmat);
void ITM_mean_norm(const double* x, const double* y, const double* h,
                   const int* n, const int* p, const int* out,
                   double* Mmat);

void ITM_pdf_norm(const double* x, const double* y, const double* hx,
                  const double* hy, const int* n, const int* p,
                  const int* out, double* Mmat);

void dlogden1(const double* x, const int* n, const int* p,
              const double* h, const int* out,
              double* den, double* dlogf);
void dlogden3(const double* x, const int* n, const int* p,
              const double* h, const int* out,
              double* den, double* rdlogf, double* dlogf);

void FM_mean(const double* x, const double* gk, const double* y,
              const double* h, const int* n, const int* p,
              const int* out, double* Mmat);

void FM_pdf(const double* x, const double* gk, const double* y,
             const double* hx, const double* hy,
             const int* n, const int* p,
             const int* out, double* Mmat);
void FM_mean_norm(const double* x, const double* y, const double* h,
                   const int* n, const int* p, const int* out,
                   double* Mmat);

void FM_pdf_norm(const double* x, const double* y, const double* hx,
                  const double* hy, const int* n, const int* p,
                  const int* out, double* Mmat);

void Fdlogden1(const double* x, const int* n, const int* p,
              const double* h, const int* out,
              double* den, double* dlogf);
void Fdlogden3(const double* x, const int* n, const int* p,
              const double* h, const int* out,
              double* den, double* rdlogf, double* dlogf);
