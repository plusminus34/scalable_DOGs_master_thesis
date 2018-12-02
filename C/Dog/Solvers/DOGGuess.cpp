#include "DOGGuess.h"

#include "GeneralizedProcrustes.h"
#include "../../Optimization/CompositeConstraints.h"

using namespace igl;

DOGGuess::DOGGuess(const Dog& dog, const bool& align_procrustes, const bool& deform_arap) : dog_init(dog), align_procrustes(align_procrustes),
															deform_arap(deform_arap) {
	Ftri = Fsqr_to_F(dog.getF());
  Vref = dog.getV();
	arapData.max_iter = 5;	
}

void DOGGuess::guess(Dog& dog, const PositionalConstraints& posConst, StitchingConstraints& stitchConst,
		EdgePointConstraints& edgePointConstraints) {
	if (align_procrustes){
		GeneralizedProcrustes genProc; genProc.solve(dog, posConst, stitchConst);	
	}
	if (deform_arap) {
		guessARAP(dog, posConst, stitchConst, edgePointConstraints);
	}
	dog.update_Vren();
}

void DOGGuess::guessARAP(Dog& dog, const PositionalConstraints& postConst, 
					StitchingConstraints& stitchConst,
					EdgePointConstraints& edgePointConstraints) {
  int vn = dog.getV().rows();
	auto b = postConst.getPositionIndices(); auto bc = postConst.getPositionVals();
	Eigen::VectorXi b_V(b.rows()/3); Eigen::MatrixXd bc_V(b_V.rows(),3);
	for (int i =0 ; i < b_V.rows(); i++) {b_V(i) = b(i);}
	vec_to_mat2(bc, bc_V);

  
  CompositeConstraints compConst({&stitchConst,&edgePointConstraints});
  //CompositeConstraints compConst({&stitchConst});
  auto x = dog.getV_vector();
  auto JacobianIJV(compConst.JacobianIJV(x));
  
  std::vector<Eigen::Triplet<double> > jacobianIJV_V(JacobianIJV.size()/3);
  
  int const_i = -1; int prev_vec_const_i = -2;
  for (int i = 0; i < JacobianIJV.size(); i++) {
    //std::cout << "JacobianIJV[i].col() = " << JacobianIJV[i].col() << std::endl;
    if (JacobianIJV[i].col() < vn) {
      // Various IJV can refer to the same constraint, and this is how we determine when to increment const_i
      if (prev_vec_const_i!= JacobianIJV[i].row()) {
        const_i++;
        prev_vec_const_i = JacobianIJV[i].row();
      }
      //std::cout << "JacobianIJV[i].row() = " << JacobianIJV[i].row() << " and const_i = " << const_i << std::endl;
      jacobianIJV_V[const_i] = Eigen::Triplet<double>(const_i,JacobianIJV[i].col(),JacobianIJV[i].value());
    }
  }
  //std::cout << "JacobianIJV_V.size() = " << jacobianIJV_V.size() << " JacobianIJV.size()/3 = " << JacobianIJV.size()/3 << std::endl;
  //std::cout << "const_num = " << const_i+1 << " compConst.getConstNum()/3 = " << compConst.getConstNum()/3 << std::endl;
  
  Eigen::SparseMatrix<double> Jacobian(compConst.getConstNum()/3, vn);
  Jacobian.setFromTriplets(jacobianIJV_V.begin(),jacobianIJV_V.end());
  
  //auto aeq(compConst.Jacobian(x));
 
  arap_precomputation_linear_equalities(Vref,Ftri,3,b_V,Jacobian,arapData);
  
  std::cout << "Jacobian = " << Jacobian << std::endl;
  std::cout << "stitchConst.Vals(x) = " << stitchConst.Vals(x) << std::endl;
  std::cout << "edgePointConstraints.Vals(x) = " << edgePointConstraints.Vals(x) << std::endl;
  Eigen::VectorXd eq_vals(compConst.Vals(x));
  Eigen::MatrixXd eq_vals_V; vec_to_mat2(eq_vals,eq_vals_V);
  std::cout << "eq_vals = " << eq_vals << std::endl;
  std::cout << "eq_vals_V = " << eq_vals_V << std::endl;
  
  arap_solve_linear_constraints(bc_V,eq_vals_V,arapData,dog.getVMutable());
  //igl::arap_precomputation(Vref,Ftri,3,b_V,arapData);
	//igl::arap_solve(bc_V,arapData,dog.getVMutable());

}

template <
  typename DerivedV,
  typename DerivedF,
  typename Derivedb>
IGL_INLINE bool DOGGuess::arap_precomputation_linear_equalities(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const int dim,
  const Eigen::PlainObjectBase<Derivedb> & b,
  const Eigen::SparseMatrix<double>& Aeq,
  ARAPData & data)
{
  using namespace std;
  using namespace Eigen;
  typedef typename DerivedV::Scalar Scalar;
  // number of vertices
  const int n = V.rows();
  data.n = n;
  std::cout << "n = " << n << std::endl;
  assert((b.size() == 0 || b.maxCoeff() < n) && "b out of bounds");
  assert((b.size() == 0 || b.minCoeff() >=0) && "b out of bounds");
  // remember b
  data.b = b;
  //assert(F.cols() == 3 && "For now only triangles");
  // dimension
  //const int dim = V.cols();
  assert((dim == 3 || dim ==2) && "dim should be 2 or 3");
  data.dim = dim;
  //assert(dim == 3 && "Only 3d supported");
  // Defaults
  data.f_ext = MatrixXd::Zero(n,data.dim);

  assert(data.dim <= V.cols() && "solve dim should be <= embedding");
  bool flat = (V.cols() - data.dim)==1;

  DerivedV plane_V;
  DerivedF plane_F;
  typedef SparseMatrix<Scalar> SparseMatrixS;
  SparseMatrixS ref_map,ref_map_dim;
  if(flat)
  {
    project_isometrically_to_plane(V,F,plane_V,plane_F,ref_map);
    repdiag(ref_map,dim,ref_map_dim);
  }
  const PlainObjectBase<DerivedV>& ref_V = (flat?plane_V:V);
  const PlainObjectBase<DerivedF>& ref_F = (flat?plane_F:F);
  SparseMatrixS L;
  cotmatrix(V,F,L);

  ARAPEnergyType eff_energy = data.energy;
  if(eff_energy == ARAP_ENERGY_TYPE_DEFAULT)
  {
    switch(F.cols())
    {
      case 3:
        if(data.dim == 3)
        {
          eff_energy = ARAP_ENERGY_TYPE_SPOKES_AND_RIMS;
        }else
        {
          eff_energy = ARAP_ENERGY_TYPE_ELEMENTS;
        }
        break;
      case 4:
        eff_energy = ARAP_ENERGY_TYPE_ELEMENTS;
        break;
      default:
        assert(false);
    }
  }


  // Get covariance scatter matrix, when applied collects the covariance
  // matrices used to fit rotations to during optimization
  covariance_scatter_matrix(ref_V,ref_F,eff_energy,data.CSM);
  if(flat)
  {
    data.CSM = (data.CSM * ref_map_dim.transpose()).eval();
  }
  assert(data.CSM.cols() == V.rows()*data.dim);

  // Get group sum scatter matrix, when applied sums all entries of the same
  // group according to G
  SparseMatrix<double> G_sum;
  if(data.G.size() == 0)
  {
    if(eff_energy == ARAP_ENERGY_TYPE_ELEMENTS)
    {
      speye(F.rows(),G_sum);
    }else
    {
      speye(n,G_sum);
    }
  }else
  {
    // groups are defined per vertex, convert to per face using mode
    if(eff_energy == ARAP_ENERGY_TYPE_ELEMENTS)
    {
      Eigen::Matrix<int,Eigen::Dynamic,1> GG;
      MatrixXi GF(F.rows(),F.cols());
      for(int j = 0;j<F.cols();j++)
      {
        Matrix<int,Eigen::Dynamic,1> GFj;
        slice(data.G,F.col(j),GFj);
        GF.col(j) = GFj;
      }
      mode<int>(GF,2,GG);
      data.G=GG;
    }
    //printf("group_sum_matrix()\n");
    group_sum_matrix(data.G,G_sum);
  }
  SparseMatrix<double> G_sum_dim;
  repdiag(G_sum,data.dim,G_sum_dim);
  assert(G_sum_dim.cols() == data.CSM.rows());
  data.CSM = (G_sum_dim * data.CSM).eval();


  arap_rhs(ref_V,ref_F,data.dim,eff_energy,data.K);
  if(flat)
  {
    data.K = (ref_map_dim * data.K).eval();
  }
  assert(data.K.rows() == data.n*data.dim);

  SparseMatrix<double> Q = (-L).eval();

  if(data.with_dynamics)
  {
    const double h = data.h;
    assert(h != 0);
    SparseMatrix<double> M;
    massmatrix(V,F,MASSMATRIX_TYPE_DEFAULT,data.M);
    const double dw = (1./data.ym)*(h*h);
    SparseMatrix<double> DQ = dw * 1./(h*h)*data.M;
    Q += DQ;
    // Dummy external forces
    data.f_ext = MatrixXd::Zero(n,data.dim);
    data.vel = MatrixXd::Zero(n,data.dim);
  }

  return min_quad_with_fixed_precompute(
    Q,b,Aeq,true,data.solver_data);
}

template <
  typename Derivedbc,
  typename DerivedU>
IGL_INLINE bool DOGGuess::arap_solve_linear_constraints(
  const Eigen::PlainObjectBase<Derivedbc> & bc,
  const Eigen::RowVectorXd& linear_const_vals,
  ARAPData & data,
  Eigen::PlainObjectBase<DerivedU> & U)
{
  using namespace Eigen;
  using namespace std;
  assert(data.b.size() == bc.rows());
  if(bc.size() > 0)
  {
    assert(bc.cols() == data.dim && "bc.cols() match data.dim");
  }
  const int n = data.n;
  int iter = 0;
  if(U.size() == 0)
  {
    // terrible initial guess.. should at least copy input mesh
#ifndef NDEBUG
    cerr<<"arap_solve: Using terrible initial guess for U. Try U = V."<<endl;
#endif
    U = MatrixXd::Zero(data.n,data.dim);
  }else
  {
    assert(U.cols() == data.dim && "U.cols() match data.dim");
  }
  // changes each arap iteration
  MatrixXd U_prev = U;
  // doesn't change for fixed with_dynamics timestep
  MatrixXd U0;
  if(data.with_dynamics)
  {
    U0 = U_prev;
  }
  while(iter < data.max_iter)
  {
    U_prev = U;
    // enforce boundary conditions exactly
    for(int bi = 0;bi<bc.rows();bi++)
    {
      U.row(data.b(bi)) = bc.row(bi);
    }

    const auto & Udim = U.replicate(data.dim,1);
    assert(U.cols() == data.dim);
    // As if U.col(2) was 0
    MatrixXd S = data.CSM * Udim;
    // THIS NORMALIZATION IS IMPORTANT TO GET SINGLE PRECISION SVD CODE TO WORK
    // CORRECTLY.
    S /= S.array().abs().maxCoeff();

    const int Rdim = data.dim;
    MatrixXd R(Rdim,data.CSM.rows());
    if(R.rows() == 2)
    {
      fit_rotations_planar(S,R);
    }else
    {
      fit_rotations(S,true,R);
//#ifdef __SSE__ // fit_rotations_SSE will convert to float if necessary
//      fit_rotations_SSE(S,R);
//#else
//      fit_rotations(S,true,R);
//#endif
    }
    //for(int k = 0;k<(data.CSM.rows()/dim);k++)
    //{
    //  R.block(0,dim*k,dim,dim) = MatrixXd::Identity(dim,dim);
    //}


    // Number of rotations: #vertices or #elements
    int num_rots = data.K.cols()/Rdim/Rdim;
    // distribute group rotations to vertices in each group
    MatrixXd eff_R;
    if(data.G.size() == 0)
    {
      // copy...
      eff_R = R;
    }else
    {
      eff_R.resize(Rdim,num_rots*Rdim);
      for(int r = 0;r<num_rots;r++)
      {
        eff_R.block(0,Rdim*r,Rdim,Rdim) =
          R.block(0,Rdim*data.G(r),Rdim,Rdim);
      }
    }

    MatrixXd Dl;
    if(data.with_dynamics)
    {
      assert(data.M.rows() == n &&
        "No mass matrix. Call arap_precomputation if changing with_dynamics");
      const double h = data.h;
      assert(h != 0);
      //Dl = 1./(h*h*h)*M*(-2.*V0 + Vm1) - fext;
      // data.vel = (V0-Vm1)/h
      // h*data.vel = (V0-Vm1)
      // -h*data.vel = -V0+Vm1)
      // -V0-h*data.vel = -2V0+Vm1
      const double dw = (1./data.ym)*(h*h);
      Dl = dw * (1./(h*h)*data.M*(-U0 - h*data.vel) - data.f_ext);
    }

    VectorXd Rcol;
    columnize(eff_R,num_rots,2,Rcol);
    VectorXd Bcol = -data.K * Rcol;
    assert(Bcol.size() == data.n*data.dim);
    for(int c = 0;c<data.dim;c++)
    {
      VectorXd Uc,Bc,bcc,Beq;
      Bc = Bcol.block(c*n,0,n,1);
      if(data.with_dynamics)
      {
        Bc += Dl.col(c);
      }
      if(bc.size()>0)
      {
        bcc = bc.col(c);
      }
      if (linear_const_vals.size()>0)
      {
        Beq = linear_const_vals.col(c);
      }
      min_quad_with_fixed_solve(
        data.solver_data,
        Bc,bcc,Beq,
        Uc);
      U.col(c) = Uc;
    }

    iter++;
  }
  if(data.with_dynamics)
  {
    // Keep track of velocity for next time
    data.vel = (U-U0)/data.h;
  }

  return true;
}