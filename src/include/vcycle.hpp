#pragma once

#include "smoothers.hpp"
#include "vector.hpp"
#include "matrices.hpp"
#include "level.hpp"

template<class T, class M, class U>
class Vcycle {
	public:
		int max_levels;
		Level<T,U,M> fine_level;
		std::vector<Levels<T,U,M>> coarse_levels;

		Vcycle(int max_levels, Level<T,U,M>& fine_level, std::vector<Levels>& mg_levels) :
			max_levels(max_levels),
			fine_level(fine_level),
			mg_levels(mg_levels)
		{}

		~Vcycle(){}
	
		void cycle(int& level) 
		{
			if (level == max_levels - 1)
			{
				coarse_levels[max_levels-1].S(coarse_levels[max_levels-1].b,coarse_level[max_levels-1].x);
				level -= 1;
				return

			}

			if (level == 0) { //on the finest level
				
				//pre smooth
				fine_level.S(fine_level.b,fine_level.x);

				//compute residual
				fine_level.A(fine_level.b,fine_level.Ax);
				fine_level.r = fine_level.b - fine_level.Ax;

				//restrict
				fine_level.R(fine_level.r,coarse_levels[0].b);

				level += 1;

				//recursive call
				this->cycle(level);

				//prolong the error to the fine level
				fine_level.P(coarse_levels[0].x,fine_level.e);

				//add it to the current solution on the fine level
				fine_level.x0 = fine_level.x + fine_level.e;

				//post smoothing
				fine_level.S(fine_level.b,fine_level.xf);

			} else {

				//pre smooth
				coarse_levels[level-1].S(coarse_levels[level-1].b,coarse_levels[level-1].x);

				//compute residual
				coarse_levels[level-1].A(coarse_levels[level-1].x,coarse_levels[level-1].Ax);

				//restrict to the next coarsest level
				coarse_levels[level-1].R(coarse_levels[level-1].r,coarse_levels[level-2].b);

				level +=1;

				//recursive call
				this->cycle(level);

				//prolong the error to the next level
				coarse_levels[level-1].P(coarse_levels[level-2].x,coarse_levels[level-1].e);

				//add it to the current solution
				coarse_levels[level-1].x0 = coarse_levels[level-1].x + coarse_levels[level-1].e;

				//post smoothing
				coarse_levels[level-1].S(coarse_levels[level-1].b,coarse_levels[level-1].xf);
			}

		}

		void operator()(int max_iters, T target)
		{

			int iter = 0;
			int level = 0;
			T rn = std::sqrt(fine_level.b.norm2());
			while (iter <= max_iters && rn >= target) {
				this->cycle(level);
				iter += 1;
				fine_level.A(fine_level.xf,fine_level.Ax);
				fine_level.r = fine_level.b - fine_level.Ax;
				rn = std::sqrt(fine_level.r.norm2());
				fine_level.x = fine_level.xf;
			}


		}



};

