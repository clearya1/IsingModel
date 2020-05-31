#include <vector>

template <class T>
class WormGrid_2_3D
{
    private:
        int 		        n_;
        std::vector<T>	    data_;

        double              beta_;
        int                 uvwz_count_;
        int                 any_2_count_;

        int index_(int x, int y, int z, int k) const
        {
            x=(x%n_+n_)%n_;
            y=(y%n_+n_)%n_;
            z=(z%n_+n_)%n_;
            return (3 * x) + (3 * n_ * y) + (3 * n_ * n_ * z) + k;
        }
    public:
        WormGrid_2_3D(int n, double beta) : n_(n), data_(n*n*n*3)
        {
            for (int i=0;i<n*n*n*3;i++) data_[i] = 0.0;
            beta_ = beta;
            uvwz_count_ = 0;
            any_2_count_ = 0;
        }

        T& operator() (int x, int y, int z, int k)
        {
            return data_[index_(x,y,z,k)];
        }
        T  operator() (int x, int y, int z, int k) const
        {
            return data_[index_(x,y,z,k)];
        }

        int size()    const { return n_; }
        
        double beta() const { return beta_; }
    
        int & uvwz_count()    { return uvwz_count_; }
    
        int & any_2_count()    { return any_2_count_; }

    
};
