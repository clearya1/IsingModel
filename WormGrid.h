#include <vector>

template <class T>
class WormGrid
{
    private:
        int 		        n_;
        std::vector<T>	    data_;

        double              beta_;
        int                 uv_count_;

        int index_(int x, int y, int z) const
        {
            x=(x%n_+n_)%n_;
            y=(y%n_+n_)%n_;
            return 2 * x + 2 * n_ * y + z;
        }
    public:
        WormGrid(int n, double beta) : n_(n), data_(n*n*2)
        {
            for (int i=0;i<n*n*2;i++) data_[i] = 0.0;
            beta_ = beta;
            uv_count_ = 0;
        }

        T& operator() (int x, int y, int z)
        {
            return data_[index_(x,y,z)];
        }
        T  operator() (int x, int y, int z) const
        {
            return data_[index_(x,y,z)];
        }

        int size()    const { return n_; }
        
        double beta() const { return beta_; }
    
        int & uv_count()    { return uv_count_; }

    
};
