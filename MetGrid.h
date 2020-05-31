#include <vector>

template <class T>
class MetGrid
{
    private:
        int 		    n_;
        std::vector<T>	data_;

        double          beta_;

        int index_(int x, int y) const
        {
            x=(x%n_+n_)%n_;
            y=(y%n_+n_)%n_;
            return x + n_ * y;
        }
    public:
        MetGrid(int n, double beta) : n_(n), data_(n*n)
        {
            for (int i=0;i<n*n;i++) data_[i] = -1.0;
            beta_ = beta;
        }

        T& operator() (int x, int y)
        {
            return data_[index_(x,y)];
        }
        T  operator() (int x, int y) const
        {
            return data_[index_(x,y)];
        }

        int size()    const { return n_; }
        
        double beta() const { return beta_; }
    
        std::vector<double> e_list;
        std::vector<double> e_sq_list;
        std::vector<double> mag_list;
        std::vector<double> mag_sq_list;
    
};
