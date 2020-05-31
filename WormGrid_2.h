#include <vector>

template <class T>
class WormGrid_2
{
    private:
        int 		        n_;
        std::vector<T>	    data_;

        double              beta_;
        int                 uvwz_count_;
        int                 any_2_count_;

        int index_(int x, int y, int z) const
        {
            x=(x%n_+n_)%n_;
            y=(y%n_+n_)%n_;
            return 2 * x + 2 * n_ * y + z;
        }
    public:
        WormGrid_2(int n, double beta) : n_(n), data_(n*n*2)
        {
            for (int i=0;i<n*n*2;i++) data_[i] = 0.0;
            beta_ = beta;
            uvwz_count_ = 0;
            any_2_count_ = 0;
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
    
        int & uvwz_count()    { return uvwz_count_; }
    
        int & any_2_count()    { return any_2_count_; }

    
};
