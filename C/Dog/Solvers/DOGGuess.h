#pragma once

#include <igl/arap.h>

#include "../Dog.h"
#include "../Objectives/StitchingConstraints.h"
#include "../../Optimization/PositionalConstraints.h"
#include "../../Optimization/EdgePointConstraints.h"

class DOGGuess {
public:
	DOGGuess(const Dog& dog, const bool& align_procrustes, const bool& deform_arap);
	void guess(Dog& dog, const PositionalConstraints& posConst, StitchingConstraints& stitchConst,
		EdgePointConstraints& edgePointConstraints);

	void update_ref(const Eigen::MatrixXd& Vref_i) {Vref = Vref_i;}
private:
	void guessARAP(Dog& dog, const PositionalConstraints& postConst, StitchingConstraints& stitchConst,
					EdgePointConstraints& edgePointConstraints);

	 template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedb>
  IGL_INLINE bool arap_precomputation_linear_equalities(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    const int dim,
    const Eigen::PlainObjectBase<Derivedb> & b,
    const Eigen::SparseMatrix<double>& Aeq,
    igl::ARAPData & data);

  template <
  typename Derivedbc,
  typename DerivedU>
IGL_INLINE bool arap_solve_linear_constraints(
  const Eigen::PlainObjectBase<Derivedbc> & bc,
  const Eigen::RowVectorXd& linear_const_vals,
  igl::ARAPData & data,
  Eigen::PlainObjectBase<DerivedU> & U);

	const bool& align_procrustes;
	const bool& deform_arap;

	const Dog& dog_init;
	Eigen::MatrixXd Vref;

	// ARAP related
	Eigen::MatrixXi Ftri; // Triangular mesh for ARAP
	igl::ARAPData arapData;
};