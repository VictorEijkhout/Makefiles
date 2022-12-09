/*  
 >
 > Stampede2 test file: This file should test that all your subroutines run.
 >
*/


/*
 * Copyright (C) 2009 - 2017 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
*/

#include <functional>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/full_matrix.h>
// IndexSet is used to set the size of each PETScWrappers::MPI::Vector:
#include <deal.II/base/index_set.h>
// PETSc appears here because SLEPc depends on this library:
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
// And then we need to actually import the interfaces for solvers that SLEPc
// provides:
#include <deal.II/lac/slepc_solver.h>
// We also need some standard C++:
#include <fstream>
#include <iostream>
//--------------for potential-------------
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_tools.h>
//------- additional includes for refinement -------
#include <deal.II/grid/grid_out.h>	// To write out refined grid files in each step
#include <deal.II/lac/constraint_matrix.h> // To handle the hanging node constraints
#include <deal.II/grid/grid_refinement.h> // to flag cells for refinement based on error indicators
#include <deal.II/numerics/error_estimator.h> //To compute the refinement indicators based on some error estimates

#include <deal.II/grid/grid_tools.h>  // To use GridTools
//------------------------------------------------

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>


#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <fstream>
#include <iostream>
#include <time.h>


// includes for gsl minimization routine
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <iomanip>                    // i/p o/p manipulator
#include <gsl/gsl_cblas.h>



namespace Step36
{
  using namespace dealii;

template <int dim>
class EigenvalueProblem
{
  public:
   EigenvalueProblem (const std::string &prm_file);
	
	double Quadrupole_Boundary_value_V_eff(const Point<dim> &p);
	double Quadrupole_Boundary_value(const Point<dim> &p);
	void run ();
	
	Vector<double> local_rho_current;
	Vector<double> local_rho_next;
	Vector<double> local_rho_previous;
	Vector<double> rho_trial;
	int scf_iteration = 0; 
	int no_electrons = 6;
	const double sigma = 0.5;
	const double *pass2rhs;
//	double Quadrupole_bv = 0.0;
	int XC_model_type=0;
	double Z_A = 6.0;

	std::vector<Point<dim>>	nuclei_locations;


	private:
	
	void move_mesh();
	void boundary_terms();
	double v_short_range(const Point<dim> &p);
	double v_external(const Point<dim> &p);
	void make_grid_and_dofs ();
	void setup_eigen_system();

	void assemble_system ();
	unsigned int solve ();
  	void output_results () const;
	void energy();
	void refine_grid();
	void E_xc_n_Hartree();
	void E_xc_SCE();

	const double E_tol = 0.0001;
	double delta_E;
	void density_mixing();

	void V_xc_LDA();
	void V_xc_SCE();

	int Nelder_Mead(unsigned int &n_b_f);
	double opt_func2 (const gsl_vector *variables_NM); 
	double v_rho_integral(std::vector<double> &std_dev_n_weights);

	std::vector<double> para_sigma_QN;
	double Quasi_Newton(std::vector<double> &parameter_sigma);
	double min_func(const gsl_vector *variables_QN);
	void min_func_grad(const gsl_vector *variables_QN, 
							 gsl_vector *grad_func);
	void min_func_n_grad(const gsl_vector *x, 
								double *f, 
								gsl_vector *grad_func);


	MPI_Comm mpi_communicator;
// no. of processes and current mpi process
	const unsigned int n_mpi_processes;
	const unsigned int this_mpi_process;

	ConditionalOStream pcout;
	
   Triangulation<dim> triangulation;
	FE_Q<dim>          fe; 
	FE_Q<dim>	fe_poisson;
   DoFHandler<dim>    dof_handler;
	DoFHandler<dim>	dof_handler_poisson;

	PETScWrappers::MPI::SparseMatrix mass_matrix_2;
	SparseMatrix<double> laplace_matrix;

   PETScWrappers::MPI::SparseMatrix		stiffness_matrix, mass_matrix;
   std::vector<PETScWrappers::MPI::Vector> eigenfunctions;
	std::vector<PETScWrappers::MPI::Vector> eig_funcs_probability;

//   std::vector<double>	eig_funcs_L2_values;
   std::vector<double>	eigenvalues;
   std::vector<double>	energy_values;
	std::vector<double>	E_Hartree;
	std::vector<double>	E_xchange_correlation;
	std::vector<double>	E_v_xc;
	std::vector<double>	V_ee_SCE; 
	std::vector<double> optimum_NM_vars;

   ParameterHandler parameters;

   ConstraintMatrix constraints;
	ConstraintMatrix constraints_rho;

	double Q_0_V_eff ;	
	Vector<double> dipole_moment_V_eff;
	FullMatrix<double> quadrupole_moment_V_eff;
	Point<dim> barycenter_V_eff;

	double Q_0 ;
	Vector<double> dipole_moment;
	FullMatrix<double>   quadrupole_moment;
	Point<dim> barycenter;

	Vector<double> psi_ground_current;
	Vector<double> psi_interpolated;
	Vector<double> epsilon_xc;
	Vector<double> epsilon_c;	
	Vector<double> epsilon_x;
	double E_xc_n_H = 0.0; 
	double E_xc = 0.0;
	double E_H = 0.0;
	double E_kin=0.0;

//--------------------------------------------------------------------------------------
//-------------- For computing the effective potential----------------------------------
// ------------- and Hartree potential
//--------------------------------------------------------------------------------------
	ConstraintMatrix constraints_V_eff;
	ConstraintMatrix constraints_V_hartree;

  	void assemble_system_V_eff ();
  	void solve_V_eff ();
	void output_results_V_eff();

	void refine_grid_apriori();
	void setup_system_V_hartree();
	void assemble_system_V_hartree();
	void solve_V_hartree();
	void output_results_V_hartree_grid (const unsigned int cycle) const;
	void output_results_V_hartree ();

	SparsityPattern 	sparsity_pattern;
	SparsityPattern sparsity_pattern_mm2; // sparsity pattern for mass matrix	
	SparseMatrix<double> system_matrix_V_eff;
	Vector<double> solution_V_eff;
	Vector<double> system_rhs_V_eff;

	PETScWrappers::MPI::SparseMatrix system_matrix_V_hartree;
	PETScWrappers::MPI::Vector solution_V_hartree;
	PETScWrappers::MPI::Vector system_rhs_V_hartree;

	Vector<double> v_xc;
   unsigned int refine_cycles = 3;  // number of refinement cycles desired.
   unsigned int ref_cycle_no = 0 ;
};




template <int dim>
EigenvalueProblem<dim>::EigenvalueProblem (const std::string &prm_file)
   :
	mpi_communicator(MPI_COMM_WORLD),
	n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator)),
	this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator)),
	pcout(std::cout, (this_mpi_process==0)),	
   fe (2),
	fe_poisson(2),	
   dof_handler (triangulation),
	dof_handler_poisson (triangulation)
{
	parameters.declare_entry ("Global mesh refinement steps", "4",
                              Patterns::Integer (0, 20),
                              "The number of times the 1-cell coarse mesh should "
                              "be refined globally for our computations.");
   parameters.declare_entry ("Number of eigenvalues/eigenfunctions", "5",
                              Patterns::Integer (0, 100),
                              "The number of eigenvalues/eigenfunctions "
                              "to be computed.");
   parameters.declare_entry ("Potential",
										"set Potential = if (x^2 + y^2 < 0.75^2, if (x*y > 0, -100, -5), 0)",
                              Patterns::Anything(),
                              "A functional description of the potential.");
 
	parameters.parse_input (prm_file);

	nuclei_locations.push_back({0.0,0.0,0.0});
}

// ------------------ Gaussian shift class for mesh motion----------------------
template <int dim>
class Gaussian_shift :public Function<dim>
{
public:
	Gaussian_shift(const Point<dim> &p);
	Point<dim> fixed_p;
	virtual void vector_value (const Point<dim> &points,
                             	Vector<double>   &values) const;
	virtual void vector_value_list  (const std::vector<Point<dim> > &points,
                                 	std::vector<Vector<double> >   &value_list) const;

};

// The constructor for this class
template<int dim>
Gaussian_shift<dim>::Gaussian_shift(const Point<dim> &p):Function<dim>(dim) 
{
	for(unsigned int i=0;i<dim;++i)
	fixed_p[i]=p[i];
}

template<int dim> 
inline
void Gaussian_shift<dim>::vector_value	(const Point<dim> &points,
                                       Vector<double>   &values) const
{
	Assert (values.size() == dim,
          ExcDimensionMismatch (values.size(), dim));
	Assert (dim >= 2, ExcNotImplemented());

	double f =1.0;
	for(unsigned int i=0;i<dim;++i)
	f *= exp(-pow(fixed_p(i)-points(i),2)/50);
		
	for(unsigned int j=0;j<dim;++j)
	values[j] = f;

}


template <int dim>
void Gaussian_shift<dim>::vector_value_list (const std::vector<Point<dim> > &points,
															std::vector<Vector<double>> &value_list) const
{
	Assert (value_list.size() == points.size(),
          ExcDimensionMismatch (value_list.size(), points.size()));
	const unsigned int n_points = points.size();
	for (unsigned int i =0; i < n_points; ++i)
		Gaussian_shift<dim>::vector_value(points[i], value_list[i]);

}



template <int dim>
void EigenvalueProblem<dim>::move_mesh ()
{

	std::vector<unsigned int> closest_vertices(nuclei_locations.size(),0);
	std::vector<bool> vertex_moved (triangulation.n_vertices(),
												false);

	for (unsigned int nuclei_no=0; nuclei_no < nuclei_locations.size(); ++nuclei_no)
	{
		const Point<dim> current_nucleus_location = nuclei_locations[nuclei_no];
		pcout<<"current nucleus location:["<<current_nucleus_location[0]<<",\t"
															<<current_nucleus_location[1]<<",\t"
															<<current_nucleus_location[2]<<"]"<<std::endl;
																					
		typename Triangulation<dim>::active_cell_iterator cell_with_nucleus = 
		GridTools::find_active_cell_around_point(triangulation,current_nucleus_location); 

		unsigned int closest_vertex = cell_with_nucleus->vertex_index(0);
		Point<dim>	closest_vertex_position = cell_with_nucleus->vertex(0);

		Point<dim> displacement_vector;

		for(unsigned int d=0; d< dim;++d)
		displacement_vector[d] = -closest_vertex_position[d] + current_nucleus_location[d];

		double smallest_distance = displacement_vector.square();
		
		for(unsigned int i =1 ; i < GeometryInfo<dim>::vertices_per_cell; ++i)
		{
			Point<dim> vertex_to_nuclei;
			for(unsigned int d=0; d< dim;++d)
			vertex_to_nuclei[d] = -cell_with_nucleus->vertex(i)[d] + current_nucleus_location[d];
			double distance = vertex_to_nuclei.square();

			if (distance < smallest_distance )
			{
				closest_vertex = cell_with_nucleus->vertex_index(i);
				closest_vertex_position = cell_with_nucleus->vertex(i);
				smallest_distance = distance;
				displacement_vector = vertex_to_nuclei;
			}


		}
		closest_vertices[nuclei_no] = closest_vertex;
		pcout<< "Closest vertex_position: \t["<< closest_vertex_position[0]<<",\t"
													<< closest_vertex_position[1]<<"\t"
													<< closest_vertex_position[2]<<"] \t closest vertex index"
													<< closest_vertex
										/*	<< closest_vertex_position[2]*/<< std::endl;



// Shifting the vertices in a gaussian manner.
		Gaussian_shift<dim> gauss_shift(closest_vertex_position);


		typename Triangulation<dim>::active_cell_iterator
					cell = triangulation.begin_active(),
					endc = triangulation.end();


		for (;cell!=endc; ++cell)
		{

			if (cell->at_boundary())
         {
            for (unsigned int face_no=0; face_no < GeometryInfo<dim>::faces_per_cell; ++ face_no)
            {
               typename Triangulation<dim>::active_face_iterator face = cell->face(face_no);
               if(face->at_boundary())
               {
                  for (unsigned int face_vertex_no =0 ;
                      face_vertex_no < GeometryInfo<dim>::vertices_per_face; ++face_vertex_no)
                  vertex_moved[face->vertex_index(face_vertex_no)] = true;
               }
            }

            for(unsigned int face_no=0; face_no < GeometryInfo<dim>::faces_per_cell; ++ face_no)
            {
               typename Triangulation<dim>::active_face_iterator face = cell->face(face_no);
               if(face->at_boundary())
               {
                  // actually needed to implement if not true condition
               }
               else
               {
                  for(unsigned int face_vertex_no =0 ;
                   face_vertex_no < GeometryInfo<dim>::vertices_per_face; ++face_vertex_no)
                  {
                     if(vertex_moved[face->vertex_index(face_vertex_no)] == false)
                     {
                        Vector<double> values(dim);
                        gauss_shift.vector_value(face->vertex(face_vertex_no), values);
                        for (unsigned int i=0; i< dim;++i)
                        face->vertex(face_vertex_no)[i] += displacement_vector[i]*values[i];

                        vertex_moved[face->vertex_index(face_vertex_no)] = true;

                     }
                  } // end for loop
               }
            }
         }	// end boundary cell loop	


			else
			{
				for(unsigned int cell_vertex_no=0; cell_vertex_no < GeometryInfo<dim>::vertices_per_cell; ++cell_vertex_no)
				{
					if (vertex_moved[cell->vertex_index(cell_vertex_no)]== false)
					{

						Vector<double> values(dim);
						gauss_shift.vector_value(cell->vertex(cell_vertex_no), values);
						for(unsigned int i=0; i< dim; ++i)
						cell->vertex(cell_vertex_no)[i] += displacement_vector[i]*values[i];


						vertex_moved[cell->vertex_index(cell_vertex_no)] = true;
					}
				}
			}
	
			for(unsigned int v=0; v< GeometryInfo<dim>::vertices_per_cell;++v)
			{
				if (cell->vertex_index(v) == closest_vertex || (cell->vertex_index(v) == closest_vertices[0] && nuclei_no ==1))
				{
					pcout<<"shifted closest vertex to position: ["
								<<cell->vertex(v)[0]<<",\t"
								<<cell->vertex(v)[1]<<",\t"
								<<cell->vertex(v)[2]<<"]"<< std::endl;
					pcout<<"cell vertex index"<<cell->vertex_index(v)<<std::endl;
				}

			}
	

		}
		for (unsigned int l=0;l<triangulation.n_vertices();++l)
		vertex_moved[l] = false;

		for (unsigned int k=0; k< nuclei_locations.size(); ++k)
		{
			if (closest_vertices[k] != 0)
			{
				vertex_moved[closest_vertices[k]] = true;	
			}
		}
		

	}
}

//==================End of mesh motion ==========================================

template <int dim>
class RightHandSide : public Function<dim>
{
public:
	RightHandSide (const double *p, std::vector<Point<dim>> &nuclei) : Function<dim>(), 
				    para_list(p), nuclei_list(nuclei) {}

	virtual double value (const Point<dim>   &p, 
      	                const unsigned int  component = 0) const;
	double value_trial_density (const Point<dim>   &p, 
   	                   const unsigned int  component = 0) const;
	double value_ion_rho (const Point<dim> &p, 
                         const unsigned int component = 0) const;
private:
	const double *para_list;
	const std::vector<Point<dim>> nuclei_list;
	int no_elec = 6;
	int no_nuclei = 1;
};

//----------------- Boundary value classes-------------------------------------------
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
	EigenvalueProblem<dim> *evp_object;
  	BoundaryValues (EigenvalueProblem<dim> *Obj_1) : Function<dim>(),
						evp_object(Obj_1) {}
  	
	virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <int dim>
class BoundaryValues_V_eff : public Function<dim>
{
public:
	EigenvalueProblem<dim> *evp_object;
	BoundaryValues_V_eff (EigenvalueProblem<dim> *Obj_1) : Function<dim>(),
						evp_object(Obj_1) {}
	
	virtual double value (const Point<dim>	&p,
	const unsigned int component=0) const;
};

// =========================================================================


template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p, 
                                  const unsigned int /*component*/) const
{
	double sigma = para_list[0];
	double return_value1 = 0.0;
	double trial_density = 0.0;
   for(unsigned int j=0; j< nuclei_list.size(); ++j)
	{  
		double rv_1 =1.0;
		double td_1 = 1.0;  
		for (unsigned int i=0; i<dim; ++i)
		{
			rv_1	*=  std::exp(-pow(p(i)-nuclei_list[j](i),2)/(std::pow(sigma,2))); 
			td_1 *=   std::exp(-pow(p(i)-nuclei_list[j](i),2)); 
		}
		return_value1 += rv_1;
		trial_density += td_1;
   }

	return (4*M_PI*trial_density)/std::pow((M_PI),1.5) - 
	(4*M_PI*return_value1)/(std::pow((M_PI*sigma*sigma),1.5));
}
  

// This function is used for the v_ion_long_range.
template <int dim>
double RightHandSide<dim>::value_ion_rho (const Point<dim> &p, 
                                  const unsigned int /*component*/) const
{
	double sigma = para_list[0];
	double return_value1 = 0.0;
   for(unsigned int j=0; j < nuclei_list.size(); ++j)
	{  
		double rv_1 =1.0;
		for (unsigned int i=0; i<dim; ++i)
		{
			rv_1	*=  std::exp(-pow(p(i)-nuclei_list[j](i),2)/(std::pow(sigma,2))); 
		}
		return_value1 += rv_1;
   }
	return (4*M_PI*return_value1)/(std::pow((M_PI*sigma*sigma),1.5));
}



template <int dim>
double RightHandSide<dim>::value_trial_density  (const Point<dim> &p, 
                                  const unsigned int /*component*/) const
{
	double trial_rho =0.0;
   for(unsigned int j=0; j< nuclei_list.size(); ++j)
	{  
		double td_1 = 1.0;  
		for (unsigned int i=0; i<dim; ++i)
		{
			td_1 *=   std::exp(-pow(p(i)-nuclei_list[j](i),2)); 
		}
		trial_rho += td_1;
   }

	return 6.0*(trial_rho)/std::pow((M_PI),1.5);  
}


// Gives the boundary values for v_H
template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
	return  evp_object->Quadrupole_Boundary_value(p);
}
 
// Gives the boundary values for v_H + v_ion_long_range or v_eff
template <int dim>
double BoundaryValues_V_eff<dim>::value (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
//	double sigma =2.0 ;

	return evp_object->Quadrupole_Boundary_value_V_eff(p);
}



//============================================== 
template <int dim>
void EigenvalueProblem<dim>::boundary_terms()
{
	Q_0 = 0.0;
	Q_0_V_eff = 0.0;
	const QGauss<dim> quadrature_formula(8);
	FEValues<dim> fe_values_rho(fe_poisson, quadrature_formula,
										update_values |
										update_quadrature_points | 
										update_JxW_values);

	const unsigned int   q_points = quadrature_formula.size();

	std::vector<double>	rho_next_cell(q_points);
	std::vector<double>	rho_previous_cell(q_points);

	dipole_moment_V_eff.reinit(3);
	quadrupole_moment_V_eff.reinit(3,3);

	dipole_moment.reinit(3);
	quadrupole_moment.reinit(3,3);
	FullMatrix<double> quadrupole_term_1(3,3);
	FullMatrix<double> quadrupole_term_1_V_eff(3,3);	
	int kronecker_mn;

	typename DoFHandler<dim>::active_cell_iterator 
	cell = dof_handler_poisson.begin_active(),
	endc = dof_handler_poisson.end();

	const RightHandSide<dim> right_hand_side (/*give value of sigma*/ pass2rhs, nuclei_locations);

	if (scf_iteration ==0)
	{  
		const MappingQ1<dim> map_for_rho;
		std::vector<Point<dim>>  support_points_rho(dof_handler_poisson.n_dofs());
//		std::map<types::global_dof_index, Point<dim>>  support_points_rho;
		DoFTools::map_dofs_to_support_points(map_for_rho,
														dof_handler_poisson, 
														support_points_rho);
		local_rho_previous.reinit(dof_handler_poisson.n_dofs());	
	
		for(unsigned int i = 0; i< dof_handler_poisson.n_dofs(); ++i)
		{
			const Point<3> sp = support_points_rho[i];
			local_rho_previous(i) = right_hand_side.value_trial_density(sp);
		}
 		if(XC_model_type==0)
		{
			for(; cell != endc; ++cell )
			{
				fe_values_rho.reinit(cell);
				fe_values_rho.get_function_values(local_rho_previous, rho_previous_cell);
				for(unsigned int i = 0; i<q_points; ++i)
				{	
					const Point<dim> position_vector = fe_values_rho.quadrature_point(i);
					Q_0 += rho_previous_cell[i] *
							fe_values_rho.JxW(i);
					dipole_moment[0] += rho_previous_cell[i] * 
											position_vector(0) *  					
										fe_values_rho.JxW(i);
					dipole_moment[1] += rho_previous_cell[i] * 
											position_vector(1) * 
											fe_values_rho.JxW(i); 
					dipole_moment[2] += rho_previous_cell[i] * 
											position_vector(2) * 
											fe_values_rho.JxW(i); 
					for (unsigned int m=0; m < 2; ++m)
						for(unsigned int n=0; n < 2; ++n)
						{
							kronecker_mn = (m == n) ? 1 : 0 ;
							quadrupole_term_1(m,n) += rho_previous_cell[i]*
												(3.0*position_vector(m)*position_vector(n) - 
												position_vector.square() * kronecker_mn) * 
												fe_values_rho.JxW(i);
						}
				}
			}	
//			pcout<<" Q_0: The total electron charge of the system for scf== 0:\t"<<Q_0<<std::endl;
			barycenter[0] = dipole_moment[0]/ Q_0;
			barycenter[1] = dipole_moment[1]/ Q_0;
			barycenter[2] = dipole_moment[2]/ Q_0;	
			for (unsigned int m=0; m < 2; ++m)
				for(unsigned int n=0; n < 2; ++n)
				{
					kronecker_mn = (m == n) ? 1 : 0 ;
					quadrupole_moment(m,n) = quadrupole_term_1(m,n) -
											3.0*Q_0*barycenter(m)*barycenter(n) + 
											Q_0*barycenter.square() * kronecker_mn; 		
				}
		}

		else
		{

		for(; cell != endc; ++cell )
		{
			fe_values_rho.reinit(cell);
			fe_values_rho.get_function_values(local_rho_previous, rho_previous_cell);
			
			for(unsigned int i = 0; i<q_points; ++i)
			{	
				const Point<dim> position_vector = fe_values_rho.quadrature_point(i);

				Q_0 += rho_previous_cell[i] *
						fe_values_rho.JxW(i);
				Q_0_V_eff += -right_hand_side.value_ion_rho(position_vector)*
						fe_values_rho.JxW(i)/(4*M_PI);
 						
				dipole_moment[0] += rho_previous_cell[i] * 
										position_vector(0) *  					
										fe_values_rho.JxW(i);
				dipole_moment[1] += rho_previous_cell[i] * 
										position_vector(1) * 
										fe_values_rho.JxW(i); 
				dipole_moment[2] += rho_previous_cell[i] * 
										position_vector(2) * 
										fe_values_rho.JxW(i); 
	
				dipole_moment_V_eff[0] += -right_hand_side.value_ion_rho(position_vector)*
										position_vector(0)*
										fe_values_rho.JxW(i)/(4*M_PI) ;
				dipole_moment_V_eff[1] += -right_hand_side.value_ion_rho(position_vector)*
										position_vector(1)*
										fe_values_rho.JxW(i)/(4*M_PI);
				dipole_moment_V_eff[2] += -right_hand_side.value_ion_rho(position_vector)*
										position_vector(2)*
										fe_values_rho.JxW(i)/(4*M_PI);

				for (unsigned int m=0; m < 2; ++m)
					for(unsigned int n=0; n < 2; ++n)
					{
						kronecker_mn = (m == n) ? 1 : 0 ;
					
						quadrupole_term_1(m,n) += rho_previous_cell[i]*
											(3*position_vector(m)*position_vector(n) - 
											position_vector.square() * kronecker_mn) * 
											fe_values_rho.JxW(i);

						quadrupole_term_1_V_eff(m,n) += -right_hand_side.value_ion_rho(position_vector)/(4*M_PI)*
											(3*position_vector(m)*position_vector(n) - 
											position_vector.square() * kronecker_mn) * 
											fe_values_rho.JxW(i);


					}

			}
		}	
	
//		pcout<<" Q_0: The total electron charge of the system for scf== 0:\t"<<Q_0<<std::endl;
		barycenter[0] = dipole_moment[0]/ Q_0;
		barycenter[1] = dipole_moment[1]/ Q_0;
		barycenter[2] = dipole_moment[2]/ Q_0;	
			            
		barycenter_V_eff[0] = dipole_moment_V_eff[0]/Q_0_V_eff;
		barycenter_V_eff[1] = dipole_moment_V_eff[1]/Q_0_V_eff;
		barycenter_V_eff[2] = dipole_moment_V_eff[2]/Q_0_V_eff;
//		pcout<<" Q_0: The total electron charge of the system for scf== 0: \t"<<Q_0_V_eff<<std::endl;

		for (unsigned int m=0; m < 2; ++m)
			for(unsigned int n=0; n < 2; ++n)
			{
				kronecker_mn = (m == n) ? 1 : 0 ;
	
				quadrupole_moment(m,n) = quadrupole_term_1(m,n) -
										3*Q_0*barycenter(m)*barycenter(n) + 
										Q_0*barycenter.square() * kronecker_mn; 		
				
				quadrupole_moment_V_eff(m,n) = quadrupole_term_1_V_eff(m,n) - 
										3*Q_0_V_eff*barycenter_V_eff(m)*barycenter_V_eff(n) + 
										Q_0_V_eff*barycenter_V_eff.square()*kronecker_mn;			
			}
	

//		pcout<<"3 Q_0: The total electron charge of the system=\t"<<Q_0<<std::endl;
		}
	}
	else
	{	
//		pcout<<"checknorm before Q_0: "<< local_rho_next.l2_norm()<<std::endl;
		for(; cell != endc; ++cell )
		{
			fe_values_rho.reinit(cell);
			fe_values_rho.get_function_values(local_rho_next, rho_next_cell);
			
			for(unsigned int i = 0; i<q_points; ++i)
			{	
				const Point<dim> position_vector = fe_values_rho.quadrature_point(i);

				Q_0 += rho_next_cell[i] *
						fe_values_rho.JxW(i);
				Q_0_V_eff += -right_hand_side.value_ion_rho(position_vector)*
						fe_values_rho.JxW(i)/(4*M_PI);

				dipole_moment[0] += rho_next_cell[i] * 
										position_vector(0) *  					
										fe_values_rho.JxW(i);

				dipole_moment[1] += rho_next_cell[i] * 
										position_vector(1) * 
										fe_values_rho.JxW(i); 
									
				dipole_moment[2] += rho_next_cell[i] * 
										position_vector(2) * 
										fe_values_rho.JxW(i); 
	
				dipole_moment_V_eff[0] += -right_hand_side.value_ion_rho(position_vector)*
										position_vector(0)*
										fe_values_rho.JxW(i)/(4*M_PI);
				dipole_moment_V_eff[1] += -right_hand_side.value_ion_rho(position_vector)*
										position_vector(1)*
										fe_values_rho.JxW(i)/(4*M_PI);
				dipole_moment_V_eff[2] += -right_hand_side.value_ion_rho(position_vector)*
										position_vector(2)*
										fe_values_rho.JxW(i)/(4*M_PI);
		
				for (unsigned int m=0; m < 2; ++m)
					for(unsigned int n=0; n < 2; ++n)
					{
						if(m=n)
							kronecker_mn = 1;
						else
							kronecker_mn = 0;
					
						quadrupole_term_1(m,n) += rho_next_cell[i]*
											(3*position_vector(m)*position_vector(n) - 
											position_vector.square() * kronecker_mn) * 
											fe_values_rho.JxW(i);
						quadrupole_term_1_V_eff (m,n) += -right_hand_side.value_ion_rho(position_vector)/(4*M_PI)*
											(3*position_vector(m)*position_vector(n) - 
											position_vector.square() * kronecker_mn) * 
											fe_values_rho.JxW(i);

					}
			}
		}	
	
		barycenter[0] = dipole_moment[0]/ Q_0;
		barycenter[1] = dipole_moment[1]/ Q_0;
		barycenter[2] = dipole_moment[2]/ Q_0;	
	
		barycenter_V_eff[0] = dipole_moment_V_eff[0]/Q_0_V_eff;
		barycenter_V_eff[1] = dipole_moment_V_eff[1]/Q_0_V_eff;
		barycenter_V_eff[2] = dipole_moment_V_eff[2]/Q_0_V_eff;
//		pcout<<"Q_0: The total electron charge of the system for new density=\t"<<Q_0
				
//<<"\tQ_0_V_eff=\t"<<Q_0_V_eff<<std::endl;
	
		for (unsigned int m=0; m < 2; ++m)
			for(unsigned int n=0; n < 2; ++n)
			{
				if(m=n)
					kronecker_mn = 1;
				else
					kronecker_mn = 0;
					
				quadrupole_moment(m,n) = quadrupole_term_1(m,n) -
										3*Q_0*barycenter(m)*barycenter(n) + 
										Q_0*barycenter.square() * kronecker_mn; 		

				quadrupole_moment_V_eff(m,n) = quadrupole_term_1_V_eff(m,n) - 
										3*Q_0_V_eff*barycenter_V_eff(m)*barycenter_V_eff(n) + 
										Q_0_V_eff*barycenter_V_eff.square()*kronecker_mn;			
			}
	}
}



template <int dim>
double EigenvalueProblem<dim>::Quadrupole_Boundary_value(const Point<dim> &p)
{
	double boundary_potential = 0.0;
	double distance=0.0;

	distance = sqrt(pow((p[0]-barycenter(0)),2) +
				pow((p[1]-barycenter(1)),2) + 
				pow((p[2]-barycenter(2)),2));

	for (unsigned int m = 0; m<2; ++m)
		for (unsigned int n = 0; n< 2; ++n)
		{
			boundary_potential += ( 0.5*(p(m)-barycenter(m))*
									quadrupole_moment(m,n)*
									(p(n)-barycenter(n)))/ pow(distance,5);
		}
	boundary_potential += Q_0 / distance;			
	return boundary_potential;
}


template <int dim>
double EigenvalueProblem<dim>::Quadrupole_Boundary_value_V_eff(const Point<dim> &p)
{

	double boundary_potential = 0.0;
	double distance_V_eff=0.0;
	double distance = 0.0;

	distance_V_eff = sqrt(pow((p[0] - barycenter_V_eff(0)),2) +
				pow((p[1]-barycenter_V_eff(1)),2) + 
				pow((p[2]-barycenter_V_eff(2)),2));
	distance = sqrt(pow((p[0]-barycenter(0)),2) +
				pow((p[1]-barycenter(1)),2) + 
				pow((p[2]-barycenter(2)),2));

	for (unsigned int m = 0; m<2; ++m)
		for (unsigned int n = 0; n< 2; ++n)
		{
			boundary_potential += ( 0.5*(p(m)-barycenter_V_eff(m))*
									quadrupole_moment_V_eff(m,n)*
									(p(n)-barycenter_V_eff(n)))/ pow(distance_V_eff,5) + 
	 								 ( 0.5*(p(m)-barycenter(m))*
									quadrupole_moment(m,n)*
									(p(n)-barycenter(n)))/ pow(distance,5);
		}

return boundary_potential;
}
//=============================================================================

// Short range potential
template<int dim>
double EigenvalueProblem<dim>::v_short_range(const Point<dim> &p)
{
	double pot = 0.0;
	for(unsigned int i=0; i< nuclei_locations.size(); ++i)
	{
		pot += -Z_A/(p - nuclei_locations[i]).norm();
	}
	return pot;
}

template<int dim>
double EigenvalueProblem<dim>::v_external(const Point<dim> &p)
{

	double pot =0.0;
	for(unsigned int i=0; i< nuclei_locations.size(); ++i)
	{
		pot += -Z_A/(p - nuclei_locations[i]).norm(); 
	}
	return pot;
}
//==============================================================================

template<int dim>
void EigenvalueProblem<dim>::setup_eigen_system()
{
	const std::vector<IndexSet> locally_owned_dofs_per_proc =
	DoFTools::locally_owned_dofs_per_subdomain(dof_handler);

	const IndexSet locally_owned_dofs = 
	locally_owned_dofs_per_proc[this_mpi_process];

	const std::vector<IndexSet> locally_owned_dofs_per_proc_poisson =
	DoFTools::locally_owned_dofs_per_subdomain(dof_handler_poisson);

	const IndexSet locally_owned_dofs_poisson = 
	locally_owned_dofs_per_proc_poisson[this_mpi_process];




	constraints.clear();	
	DoFTools::make_hanging_node_constraints(dof_handler, constraints);	
	VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  	constraints.close ();

	constraints_rho.clear();
	DoFTools::make_hanging_node_constraints(dof_handler_poisson,constraints_rho);
	VectorTools::interpolate_boundary_values(dof_handler_poisson,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints_rho);

//	DoFTools::make_zero_boundary_constraints (dof_handler_poisson, constraints_rho);
	constraints_rho.close();

	DynamicSparsityPattern dsp_evp(locally_owned_dofs);
	DoFTools::make_sparsity_pattern (dof_handler,
												dsp_evp,
												constraints,
												false);
	std::vector<dealii::types::global_dof_index> n_locally_owned_dofs(
   	 n_mpi_processes);
	for (unsigned int i = 0; i < n_mpi_processes; i++)
   	n_locally_owned_dofs[i] = locally_owned_dofs_per_proc[i].n_elements();

  	SparsityTools::distribute_sparsity_pattern(dsp_evp,
                                              n_locally_owned_dofs,
                                              mpi_communicator,
                                              locally_owned_dofs);

    stiffness_matrix.reinit (locally_owned_dofs,
                             locally_owned_dofs,
                             dsp_evp,
										mpi_communicator);
    mass_matrix.reinit (locally_owned_dofs,
                        locally_owned_dofs,
                        dsp_evp,
								mpi_communicator);

	DynamicSparsityPattern dsp_mass(locally_owned_dofs);
	DoFTools::make_sparsity_pattern (dof_handler,
												dsp_mass,
												constraints,
												true);

	mass_matrix_2.reinit(locally_owned_dofs,
								locally_owned_dofs,
								dsp_mass,
								mpi_communicator);

   eigenfunctions.resize (parameters.get_integer ("Number of eigenvalues/eigenfunctions"));
   for (unsigned int i=0; i<eigenfunctions.size (); ++i)
   {
		eigenfunctions[i].reinit (locally_owned_dofs, mpi_communicator);
		eigenfunctions[i].compress(VectorOperation::insert);
	}

   eigenvalues.resize (eigenfunctions.size ());

   eig_funcs_probability.resize (no_electrons/2.0);

   for (unsigned int i=0; i<eig_funcs_probability.size (); ++i)
     eig_funcs_probability[i].reinit (locally_owned_dofs_poisson, 
												  mpi_communicator);

}


template <int dim>
void EigenvalueProblem<dim>::make_grid_and_dofs ()
{
	 const std::vector<IndexSet> locally_owned_dofs_per_proc =
   DoFTools::locally_owned_dofs_per_subdomain(dof_handler);

   const IndexSet locally_owned_dofs =
   locally_owned_dofs_per_proc[this_mpi_process];


	if (scf_iteration==0)
	{
		
	
		constraints.clear();	

		DoFTools::make_hanging_node_constraints(dof_handler, constraints);	
   	DoFTools::make_zero_boundary_constraints (dof_handler, constraints);
  		constraints.close ();
		
		constraints_rho.clear();
      DoFTools::make_hanging_node_constraints(dof_handler_poisson, constraints_rho);
      VectorTools::interpolate_boundary_values(dof_handler_poisson,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints_rho);
      constraints_rho.close();

	}
	else
	{
		constraints.clear();	

		DoFTools::make_hanging_node_constraints(dof_handler, constraints);	
	   DoFTools::make_zero_boundary_constraints (dof_handler, constraints);
	   constraints.close ();

		constraints_rho.clear();
      DoFTools::make_hanging_node_constraints(dof_handler_poisson, constraints_rho);
      VectorTools::interpolate_boundary_values(dof_handler_poisson,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints_rho);
      constraints_rho.close();

		solution_V_eff.reinit (dof_handler_poisson.n_dofs());
		system_rhs_V_eff.reinit (dof_handler_poisson.n_dofs());

 		constraints_V_eff.clear();

		DoFTools::make_hanging_node_constraints (dof_handler_poisson, constraints_V_eff);
		boundary_terms();
		VectorTools::interpolate_boundary_values(dof_handler_poisson,
															0,
															BoundaryValues_V_eff<dim>(this),
															constraints_V_eff);
		constraints_V_eff.close();
	
		DynamicSparsityPattern dsp_V_eff(dof_handler_poisson.n_dofs());
		DoFTools::make_sparsity_pattern (dof_handler_poisson,
													dsp_V_eff,
													constraints_V_eff, 
													false); // the false is for 'keep constrained degrees of freedom'	
		sparsity_pattern.copy_from(dsp_V_eff);
		system_matrix_V_eff.reinit(sparsity_pattern);


		constraints_V_hartree.clear();
		DoFTools::make_hanging_node_constraints(dof_handler_poisson, constraints_V_hartree);
		VectorTools::interpolate_boundary_values(dof_handler_poisson,  
															0,
															BoundaryValues<dim>(this),
															constraints_V_hartree);
		constraints_V_hartree.close();
      DynamicSparsityPattern dsp(dof_handler_poisson.n_dofs(), dof_handler_poisson.n_dofs());

		DoFTools::make_sparsity_pattern (dof_handler_poisson,
													dsp,
													constraints_V_hartree,
													true );

    const std::vector<IndexSet> locally_owned_dofs_per_proc =
      DoFTools::locally_owned_dofs_per_subdomain(dof_handler_poisson);

      const IndexSet locally_owned_dofs =
      locally_owned_dofs_per_proc[this_mpi_process];

      system_matrix_V_hartree.reinit(locally_owned_dofs,
                           locally_owned_dofs,
                           dsp,
                           mpi_communicator);

      solution_V_hartree.reinit(locally_owned_dofs, mpi_communicator);
      system_rhs_V_hartree.reinit(locally_owned_dofs, mpi_communicator);



// ------------- End V_hartree and the scf loop ---------------------
	}

	DynamicSparsityPattern dsp_evp(locally_owned_dofs);
   DoFTools::make_sparsity_pattern (dof_handler,
                                    dsp_evp,
                                    constraints,
                                    false);
   std::vector<dealii::types::global_dof_index> n_locally_owned_dofs(
       n_mpi_processes);
   for (unsigned int i = 0; i < n_mpi_processes; i++)
      n_locally_owned_dofs[i] = locally_owned_dofs_per_proc[i].n_elements();

   SparsityTools::distribute_sparsity_pattern(dsp_evp,
                                              n_locally_owned_dofs,
                                              mpi_communicator,
                                              locally_owned_dofs);

    stiffness_matrix.reinit (locally_owned_dofs,
                             locally_owned_dofs,
                             dsp_evp,
                              mpi_communicator);
    mass_matrix.reinit (locally_owned_dofs,
                        locally_owned_dofs,
                        dsp_evp,
                        mpi_communicator);

    eigenfunctions.resize (parameters.get_integer ("Number of eigenvalues/eigenfunctions"));
    for (unsigned int i=0; i<eigenfunctions.size (); ++i)
      eigenfunctions[i].reinit (locally_owned_dofs, mpi_communicator);

    eigenvalues.resize (eigenfunctions.size ());
}



template <int dim>
void EigenvalueProblem<dim>::output_results_V_eff ()
{
	if(this_mpi_process==0)
	{
		DataOut<dim> data_out_V_eff;
		data_out_V_eff.attach_dof_handler(dof_handler);
		data_out_V_eff.add_data_vector(solution_V_eff, "solution_V_eff");
		data_out_V_eff.build_patches();

		std::ofstream output_V_eff(dim == 2 ? "solution_V_kantorovich.vtk":
												  "solution_V_kantorovich.vtk");
		

		data_out_V_eff.write_vtk (output_V_eff);
	}
}


template <int dim>
void EigenvalueProblem<dim>::setup_system_V_hartree()
{
	GridTools::partition_triangulation(n_mpi_processes, triangulation);
	
	dof_handler_poisson.distribute_dofs(fe_poisson);
	DoFRenumbering::subdomain_wise (dof_handler_poisson);	

	dof_handler.distribute_dofs (fe);
	DoFRenumbering::subdomain_wise (dof_handler);	


	constraints_V_hartree.clear();
	DoFTools::make_hanging_node_constraints(dof_handler_poisson, constraints_V_hartree);
	boundary_terms();
	VectorTools::interpolate_boundary_values(dof_handler_poisson,  
															0,
															BoundaryValues<dim>(this),
															constraints_V_hartree);
	constraints_V_hartree.close();   
}


template <int dim>
void EigenvalueProblem<dim>::assemble_system_V_hartree()
{
	const QGauss<dim>  quadrature_formula_V_hartree(8);
  	const RightHandSide<dim> right_hand_side(/*give value of sigma*/ pass2rhs,nuclei_locations);

  	FEValues<dim> fe_values_V_hartree (fe_poisson,
											quadrature_formula_V_hartree,
                           		update_values    |
											update_gradients |
            				         update_quadrature_points | 
											update_JxW_values);
  
	
	const unsigned int   dofs_per_cell_V_hartree = fe_poisson.dofs_per_cell;
  	const unsigned int   n_q_points_V_hartree    = quadrature_formula_V_hartree.size();

  	FullMatrix<double>   cell_matrix_V_hartree (dofs_per_cell_V_hartree, dofs_per_cell_V_hartree);
  	Vector<double>       cell_rhs_V_hartree (dofs_per_cell_V_hartree);

	std::vector<double>	rho_next_cell(n_q_points_V_hartree);
	
  	std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell_V_hartree);

	typename DoFHandler<dim>::active_cell_iterator
  	cell_poisson = dof_handler_poisson.begin_active(),
  	endc_poisson = dof_handler_poisson.end();
	
	if(scf_iteration == 0)
	{	
		pcout<<"scf_iteration =0, assemble loop"<<std::endl;
		for (; cell_poisson!=endc_poisson; ++cell_poisson)
  		{
  			fe_values_V_hartree.reinit (cell_poisson);
			cell_matrix_V_hartree = 0;
			cell_rhs_V_hartree = 0;
   
   	   for (unsigned int q_index=0; q_index < n_q_points_V_hartree; ++q_index)
       		for (unsigned int i=0; i<dofs_per_cell_V_hartree; ++i)
          	{
            	for (unsigned int j=0; j < dofs_per_cell_V_hartree; ++j)
	      		cell_matrix_V_hartree(i,j) += (fe_values_V_hartree.shape_grad (i, q_index) *
                                   fe_values_V_hartree.shape_grad (j, q_index) *
                                   fe_values_V_hartree.JxW (q_index));
// Change value_ion_rho to value later
					cell_rhs_V_hartree(i) += 4.0*M_PI* (fe_values_V_hartree.shape_value (i, q_index) *
                            right_hand_side.value_trial_density (fe_values_V_hartree.quadrature_point (q_index)) *
                            fe_values_V_hartree.JxW (q_index));
      	    }


     		cell_poisson->get_dof_indices (local_dof_indices);
			constraints_V_hartree.distribute_local_to_global (cell_matrix_V_hartree,
															 cell_rhs_V_hartree,
															 local_dof_indices,
															 system_matrix_V_hartree,
															 system_rhs_V_hartree);  

		
//			++cell_eigen;	

		}	//end of for loop
	}
	else
	{
		for (; cell_poisson!=endc_poisson; ++cell_poisson)
  		{
  			fe_values_V_hartree.reinit (cell_poisson);
			
			cell_matrix_V_hartree = 0;
			cell_rhs_V_hartree = 0;

      	fe_values_V_hartree.get_function_values(local_rho_next, rho_next_cell);

			for (unsigned int q_index = 0; q_index < n_q_points_V_hartree; ++q_index)
        		for (unsigned int i=0; i<dofs_per_cell_V_hartree; ++i)
          	{
					for (unsigned int j=0; j < dofs_per_cell_V_hartree; ++j)
              		cell_matrix_V_hartree(i,j) += (fe_values_V_hartree.shape_grad (i, q_index) *
            			                       fe_values_V_hartree.shape_grad (j, q_index) *
				                                fe_values_V_hartree.JxW (q_index));
								
            		cell_rhs_V_hartree(i) += (fe_values_V_hartree.shape_value (i, q_index) *
                            				(4.0*M_PI*rho_next_cell[q_index]) *
                            				fe_values_V_hartree.JxW (q_index));
          	}

      	cell_poisson->get_dof_indices (local_dof_indices);
			constraints_V_hartree.distribute_local_to_global (cell_matrix_V_hartree,
															 cell_rhs_V_hartree,
															 local_dof_indices,
															 system_matrix_V_hartree,
															 system_rhs_V_hartree);  
		}
	}
}


template <int dim>
void EigenvalueProblem<dim>::solve_V_hartree ()
{
   SolverControl      solver_control (solution_V_hartree.size(),1e-12);
   PETScWrappers::SolverCG cg(solver_control, mpi_communicator);

   PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix_V_hartree);

   cg.solve (system_matrix_V_hartree, solution_V_hartree, system_rhs_V_hartree,
                preconditioner);

// distributing the constraints is tricky because of the wya the nodes are distributed
   Vector<double> temp_localized_solution_V_hartree(solution_V_hartree);
   constraints_V_hartree.distribute (temp_localized_solution_V_hartree);
   solution_V_hartree= temp_localized_solution_V_hartree;

   return solver_control.last_step();

}

template <int dim>
void EigenvalueProblem<dim>::refine_grid_apriori()
{ 
	Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

	if (ref_cycle_no < 2)
   {
      KellyErrorEstimator<dim>::estimate (dof_handler_poisson,
                                      QGauss<dim-1>(3),
                                      typename FunctionMap<dim>::type(),
                                      local_rho_previous,
                                      estimated_error_per_cell);
   }
   else
   {
      KellyErrorEstimator<dim>::estimate (dof_handler_poisson,
                                      QGauss<dim-1>(3),
                                      typename FunctionMap<dim>::type(),
                                      local_rho_current,
                                      estimated_error_per_cell);
   }
	double refine_fraction = 0.0;
   if (ref_cycle_no==0)
   {
      refine_fraction = 0.4;
   }
   else if(ref_cycle_no==1)
   {
      refine_fraction = 0.2;
   }
   else if(ref_cycle_no==2)
   {
      refine_fraction = 0.1;
   }
   else
   {
      refine_fraction = 0.08;
   }
   GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                 estimated_error_per_cell,
                                                 refine_fraction, 0.02);


	triangulation.execute_coarsening_and_refinement();
}



template <int dim>
void EigenvalueProblem<dim>::output_results_V_hartree_grid (const unsigned int cycle) const
{
  	Assert (cycle < 10, ExcNotImplemented());
	std::string filename = "grid-";
  	filename += ('0' + cycle);
  	filename += ".vtk";
  	std::ofstream output (filename.c_str());
  	GridOut grid_out;
  	grid_out.write_vtk (triangulation, output);
}


template <int dim>
void EigenvalueProblem<dim>::output_results_V_hartree ()
{
	DataOut<dim> data_out_V_hartree;
	data_out_V_hartree.attach_dof_handler(dof_handler_poisson);
	data_out_V_hartree.add_data_vector(solution_V_hartree, "solution_V_hartree");
	data_out_V_hartree.build_patches();

	std::ofstream output_V_hartree(dim == 2 ? "solution_V_hartree-2d.vtk":
											  "solution_V_hartree-3d.vtk");
	

	data_out_V_hartree.write_vtk (output_V_hartree);
}



//=================================================================================
//====================== End Of LDA ================================

template <int dim>
class Potential_function : public Function<dim>
{
public:
	const std::vector<double> *var_list;
	double N_e;
	std::vector<Point<dim>> R_I;

	Potential_function( const std::vector<double> *p,
								const double &n_particles, 
								std::vector<Point<dim>> &nuclei_locations
							) : Function<dim>(), 
				var_list(p),
				N_e(n_particles),
				R_I(nuclei_locations)
 				{ };


	const std::vector<double> &var_ref = *var_list;

  	virtual double value (const Point<dim>   &pt, 
                       const unsigned int  component = 0) const;
	double v_kantorovich (const Point<dim> &pt) const ;	
	void  v_kantorovich_grad (const Point<dim> &pt, Vector<double> &values);
};




template<int dim>
double Potential_function<dim>::value (const Point<dim> &pt, const unsigned int) const
{
	return 1.0;
}

template<int dim>
double Potential_function<dim>::v_kantorovich(const Point<dim> &pt) const
{
	int basis_size = var_list->size();
	double v_k = 0.0;

	double dr_weights =0.0;
	for (int i =0; i< 4; ++i)
	dr_weights +=  var_ref[5*i]*var_ref[5*i];// +var_ref[2]*var_ref[2];

	const Point<3> origin = {0.00000001,0,0};
	

	for (int j=0; j < 4 ; ++j)
	{
		const Point<3> var_pt = {var_ref[5*j+2],var_ref[5*j+3],var_ref[5*j+4]};
		if(pt == var_pt)
		{
			v_k += pow(var_ref[5*j],2)/dr_weights /(origin).norm() *
					 erf((origin).norm()/(sqrt(2.0)*var_ref[5*j+1])) ;
		}
		else
		{	
			v_k +=  pow(var_ref[5*j],2)/dr_weights /(pt-var_pt).norm() *
					 erf((pt-var_pt).norm()/(sqrt(2.0)*var_ref[5*j+1]));
		}
	}

	return 5.0*v_k;
}



//__________________________________________
template<int dim>
void Potential_function<dim>::v_kantorovich_grad(const Point<dim> &pt, Vector<double> &values )
{
	int basis_size = var_list->size();
	
	double dr_weights =0.0;
	for (int i=0; i < 4; ++i)
		dr_weights += var_ref[5*i]*var_ref[5*i];
	

	const Point<3> origin = {0.00000001,0,0};

	double term1 = 0.0;

// basis := [weight, sigma, Rx, Ry , Rz]

	for (int j=0; j < 4 ; ++j)
	{
		term1 = 0.0;
		const Point<3> var_pt = {var_ref[5*j+2],var_ref[5*j+3],var_ref[5*j+4]};

		if(pt==var_pt)
		{
			term1 += pow(var_ref[5*j],2)/dr_weights*(1.0/ pow((origin).norm(),3) * 
          	      erf((origin).norm()/(sqrt(2.0)*var_ref[5*j+1])) - 
		        	   2.0/sqrt(M_PI*2.0)/var_ref[5*j+1] /
						pow((origin).norm(),2) * 
						exp(- pow((origin).norm(),2)/2.0/pow(var_ref[5*j+1],2)));
		}
		else
		{
			term1 += pow(var_ref[5*j],2)/dr_weights*(1.0/ pow((pt-var_pt).norm(),3) * 
          	       erf((pt-var_pt).norm()/(sqrt(2.0)*var_ref[5*j+1])) - 
              	    2.0/sqrt(M_PI*2.0)/var_ref[5*j+1] /
						 pow((pt-var_pt).norm(),2.0) * 
	                exp(- pow((pt-var_pt).norm(),2)/2.0/pow(var_ref[5*j+1],2)));
		}
		
		values(0) += (pt[0] - var_pt[0] )*term1;// + (pt[0]-right_pt[0])*term2;
		values(1) += (pt[1] - var_pt[1] )*term1;// + (pt[1]-right_pt[1])*term2;
		values(2) += (pt[2] - var_pt[2] )*term1;// + (pt[2]-right_pt[2])*term2;
	}


	values(0) *=  5.0;
	values(1) *=  5.0;
 	values(2) *=  5.0;
}
//________________________________________



// ==== Kantorovich potential: Strictly Correlated Electron limit of DFT============
template<int dim>
void EigenvalueProblem<dim>::V_xc_SCE()
{
	unsigned int n_basis_funcs = 5;	
	Nelder_Mead(n_basis_funcs);
}


// Wrapper for gsl_multimin_function

class gsl_function_pp : public gsl_multimin_function
{
	public:
	gsl_function_pp(std::function<double(const gsl_vector *)> const& func): _func(func)
	{
		f = &gsl_function_pp::invoke;
		params = this;
	}
	private:
		std::function<double(const gsl_vector *)> _func;
		static double invoke(const gsl_vector* x, void *params)
		{
			return static_cast<gsl_function_pp*>(params)->_func(x);
		}

};

// Wrapper for gsl_multimin_function_fdf
class gsl_function_pp2 :public gsl_multimin_function_fdf
{
	public:
	gsl_function_pp2 (std::function<double(const gsl_vector *)> const& func_QN,
							std::function<void(const gsl_vector *,gsl_vector* )>const& dfunc,
							std::function<void(const gsl_vector*,double *, gsl_vector*)>const& fdfunc):
							_func_QN(func_QN), _dfunc(dfunc), _fdfunc(fdfunc)
	{
		f = &gsl_function_pp2::invoke_f;
		df= &gsl_function_pp2::invoke_df;
		fdf=&gsl_function_pp2::invoke_fdf;
		params = this;
	}
	private:
		std::function<double(const gsl_vector *)> _func_QN;
		static double invoke_f(const gsl_vector* x, void *params)
		{
			return static_cast<gsl_function_pp2*>(params)->_func_QN(x);
		}

		std::function<void(const gsl_vector *,gsl_vector* )> _dfunc;
		static void invoke_df(const gsl_vector *x,void *params, gsl_vector *g)
		{
			return static_cast<gsl_function_pp2*>(params)->_dfunc(x,g);
		}


		std::function<void(const gsl_vector*,double *, gsl_vector*)> _fdfunc;
		static void invoke_fdf(const gsl_vector *x,void *params, double *f, gsl_vector *g)
		{
			return static_cast<gsl_function_pp2*>(params)->_fdfunc(x,f,g);
		}
};


template<int dim>
int EigenvalueProblem<dim>::Nelder_Mead(unsigned int &n_b_f)
{	
//
	size_t iter =0;
	int status=0;
	double size;

	size_t n = n_b_f;     


	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;		// optimization type
	gsl_multimin_fminimizer *s; 		// pointer to the optimzation type for an n dimension function
	s  = gsl_multimin_fminimizer_alloc( T, n);
	
	gsl_vector *x, *step_size;

	x = gsl_vector_alloc(n);

	step_size = gsl_vector_alloc(n);

// 1 basis := [weight, sigma, Rx, Ry, Rz]
   if (scf_iteration==0 || scf_iteration==1)
   {
		for (int i=0; i<1; ++i)
		{
	      gsl_vector_set(x,5*i,1);
			gsl_vector_set(x,5*i+1,1.0+0.5*i);
			gsl_vector_set(x,5*i+2,0.0);
			gsl_vector_set(x,5*i+3,0.0);
			gsl_vector_set(x,5*i+4,0.0);


      	gsl_vector_set(step_size,5*i,0.5);
      	gsl_vector_set(step_size,5*i+1,0.1);
      	gsl_vector_set(step_size,5*i+2,0.1);
      	gsl_vector_set(step_size,5*i+3,0.1);
      	gsl_vector_set(step_size,5*i+4,0.1);

		}
	}
   else
   {
		for(unsigned int i=0; i< n_b_f;++i)
		{
	      gsl_vector_set(x,i,optimum_NM_vars[i]);
	      gsl_vector_set(step_size,i,0.1);
		}
   }

	gsl_function_pp Fp(std::bind(&Step36::EigenvalueProblem<dim>::opt_func2, &(*this),std::placeholders::_1));
	gsl_multimin_function *V_ee = static_cast<gsl_multimin_function*>(&Fp) ;		

	V_ee->n = n_b_f;
	
	gsl_multimin_fminimizer_set(s, V_ee, x, step_size);
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if(status)
		{
			pcout<<"**********Error in Nelder_Mead*********"<<std::endl;
			break;
		}

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size (size,1e-8);

		if(status == GSL_SUCCESS)
		{
			pcout<<"Converged to minimum: "<<std::endl;
		
		}
   }
   while (status == GSL_CONTINUE && iter< 10);

   pcout<<"GSL, Nelder Mead, Quasi Newton routine used successfully_!!\n"
            <<"\t iteration no.:"<<iter
            <<"\t The optimum value of Vee_SCE funcitonal __=__"<< s->fval
           	<<"\t size of simplex = \t"<<size
				<<std::endl;

   pcout<<"Optimum parameters:"
				<<"\n__w1 =\t "<<gsl_vector_get(s->x,0)
				<<"__\u03C31 ="<< gsl_vector_get(s->x,1)
            <<std::endl;
		

	optimum_NM_vars.clear();
	for (unsigned int j=0; j <  n_b_f; ++j)
	{
		pcout<<"Optimum parameters: "<<
					gsl_vector_get(s->x,j)<<
					std::endl;
		optimum_NM_vars.push_back(gsl_vector_get(s->x,j));
	}		
	V_ee_SCE.push_back(s->fval);

	gsl_vector_free(x);
	gsl_vector_free(step_size);
	gsl_multimin_fminimizer_free(s);

	return status;
}
//_________________________________________ .




template<int dim>
double EigenvalueProblem<dim>::opt_func2(const gsl_vector *variables_NM)
{
	std::vector<double> vars_NM;//[variables_NM->size];
	for (unsigned int i=0; i < variables_NM->size ; ++i)
		vars_NM.push_back( gsl_vector_get(variables_NM,i));

	double v_functional = 0.0;
	v_functional = -1.0*(v_rho_integral(vars_NM) + Quasi_Newton(vars_NM));
	return v_functional;
}



template <int dim>
double EigenvalueProblem<dim>::v_rho_integral(std::vector<double> &std_dev_n_weights)
{
	Potential_function<dim> v_kant_NM (&std_dev_n_weights, no_electrons, nuclei_locations);
	QGauss<dim> quadrature_formula_kant(8);
	FEValues<dim> fe_values_kant(fe_poisson, quadrature_formula_kant,
										update_values |
										update_quadrature_points |
										update_JxW_values);	

	const unsigned int n_quad_points = quadrature_formula_kant.size();

	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler_poisson.begin_active(),
	end_cell = dof_handler_poisson.end();
	
	Point<dim>	eval_point;
	double integral_ans = 0.0;

	std::vector<double> rho_next_cell(n_quad_points);
	std::vector<double> rho_previous_cell(n_quad_points);

	if (scf_iteration==0)
	{
		for(;cell != end_cell; cell++)
		{
			fe_values_kant.reinit(cell);
			fe_values_kant.get_function_values(local_rho_previous, rho_previous_cell);
			for(unsigned int i = 0; i < n_quad_points; ++i)
			{
				eval_point =  fe_values_kant.quadrature_point(i);		//quadrature_point gives the position of the i-th quadrature point in real space.
				integral_ans += v_kant_NM.v_kantorovich(eval_point) *
								rho_previous_cell[i] *
								fe_values_kant.JxW(i);
			}
		}
	}
	else
	{
		for(;cell != end_cell; cell++)
		{
			fe_values_kant.reinit(cell);
			fe_values_kant.get_function_values(local_rho_next, rho_next_cell);
			for(unsigned int i = 0; i < n_quad_points; ++i)
			{
				eval_point =  fe_values_kant.quadrature_point(i);
				integral_ans += v_kant_NM.v_kantorovich(eval_point) *
									rho_next_cell[i] *
									fe_values_kant.JxW(i);
			}
		}
	}

	return integral_ans;
}
//_____________________



template<int dim>
double EigenvalueProblem<dim>::Quasi_Newton(std::vector<double> &para_list)
{
	size_t iter = 0;
	int status = 0;

	const gsl_multimin_fdfminimizer_type *T; // minimizer of type T
	gsl_multimin_fdfminimizer *s;            // This function is used to initialize a multidimensional minimizer.


	para_sigma_QN.clear();
	para_sigma_QN = para_list;

	gsl_function_pp2 FDFp(std::bind(&Step36::EigenvalueProblem<dim>::min_func, &(*this), std::placeholders::_1),
				std::bind(&Step36::EigenvalueProblem<dim>::min_func_grad, &(*this), std::placeholders::_1,std::placeholders::_2),
				std::bind(&Step36::EigenvalueProblem<dim>::min_func_n_grad, &(*this), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
		
	gsl_multimin_function_fdf *g_v = static_cast<gsl_multimin_function_fdf*>(&FDFp);
        
	g_v->n  = 3*no_electrons;
 
	T = gsl_multimin_fdfminimizer_vector_bfgs2;
	s = gsl_multimin_fdfminimizer_alloc (T,3*no_electrons);

	gsl_vector *x = gsl_vector_alloc (3*no_electrons);
	
	for(int i=0;i<no_electrons ; ++i)
	{
	   gsl_vector_set (x , 3*i   ,0);
	   gsl_vector_set (x , 3*i+1 ,2.0*i - 5.0);
	   gsl_vector_set (x , 3*i+2 ,2.0*i - 5.0);
	}


	gsl_multimin_fdfminimizer_set (s, g_v, x, 0.1,1e-3);
	do
	{
 		iter++;


		status = gsl_multimin_fdfminimizer_iterate (s);

		if (status)
		{ 
		  break;	// This should happen only when an error occurs 
		}
		status = gsl_multimin_test_gradient (s->gradient, 1e-9);	// set tolerance
		
		if (status == GSL_SUCCESS)
		{ 

		}
	}
	while (status == GSL_CONTINUE && iter < 10);

	double Optimum_gv = s->f;
	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);
	return Optimum_gv;		
}


template<int dim>
double EigenvalueProblem<dim>::min_func(const gsl_vector *variables_QN)
{
	unsigned int N = 6;
	double r[3*N] = {0.0};
	Potential_function<dim> v_kant_QN(&para_sigma_QN, no_electrons, nuclei_locations);

	double return1=0.0;
	double return2=0.0;

	for(int i=0;i < no_electrons;++i )
	{	r[3*i]   = gsl_vector_get(variables_QN,3*i);
		r[3*i+1] = gsl_vector_get(variables_QN,3*i+1);
		r[3*i+2] = gsl_vector_get(variables_QN,3*i+2);
  	} 

	for(unsigned int i=0; i<N;++i)
   {
		for(unsigned int j = 0; j<N; ++j)
	 	{
      	if (j>i)
	      {
	      	return1 += 1.0/(sqrt(pow((r[3*i]-r[3*j]),2)+pow((r[3*i+1]-r[3*j+1]),2)
                       +pow((r[3*i+2]-r[3*j+2]),2)));
   	   }
	   	else
	     	{
	      	return1 += 0;
	     	}
	 	}
  
		const Point<3> p(r[3*i],r[3*i+1],r[3*i+2]);
		return2 = return2 +  v_kant_QN.v_kantorovich(p);

	}	        
	return return1-return2;  
}



template<int dim>
void EigenvalueProblem<dim>::min_func_grad(const gsl_vector *variables_QN, gsl_vector *grad_func)
{
   double r[3*no_electrons];
	Potential_function<dim> v_kant_QN_grad(&para_sigma_QN, no_electrons, nuclei_locations);
   for(int i=0;i < no_electrons;++i )
   {
		r[3*i]   = gsl_vector_get(variables_QN,3*i);
		r[3*i+1] = gsl_vector_get(variables_QN,3*i+1);
	  	r[3*i+2] = gsl_vector_get(variables_QN,3*i+2);
  } 
   double return3[no_electrons]={0};
   double return4[no_electrons]={0};
   double return5[no_electrons]={0};
	Vector<double> v_gradient(3);

	for(int i=0; i<no_electrons; ++i)
	{   
		for(int j=0; j<no_electrons; ++j)
	 	{
      	if (j!=i) 
	     	{
				//The x component: of the electron-electron repulsion potential
     
	      	return3[i]+= -1.0*(r[3*i]-r[3*j])/pow((sqrt(pow(r[3*i]-r[3*j],2)+pow(r[3*i+1]-r[3*j+1],2)
                      +pow(r[3*i+2]-r[3*j+2],2))),3);             
				//y-component:     
		      return4[i]+= -1.0*(r[3*i+1]-r[3*j+1])/pow((sqrt(pow(r[3*i]-r[3*j],2)+pow(r[3*i+1]-r[3*j+1],2)
                      +pow(r[3*i+2]-r[3*j+2],2))),3);     
				//z-component:   
	   	   return5[i]+= -1.0*(r[3*i+2]-r[3*j+2])/pow((sqrt(pow(r[3*i]-r[3*j],2)+pow(r[3*i+1]-r[3*j+1],2)
                      +pow(r[3*i+2]-r[3*j+2],2))),3);     
	
			}
			else
			{
     		 	return3[i] += 0;
      		return4[i] += 0;
       		return5[i] += 0;
	      }
		}

		const Point<3> p(r[3*i],r[3*i+1],r[3*i+2]);
		v_gradient.reinit(3);
		v_kant_QN_grad.v_kantorovich_grad(p,v_gradient);

		return3[i] += v_gradient(0);

		return4[i] += v_gradient(1);

		return5[i] += v_gradient(2);

      gsl_vector_set(grad_func,3*i,return3[i]);	
      gsl_vector_set(grad_func,3*i+1,return4[i]);	
      gsl_vector_set(grad_func,3*i+2,return5[i]);
	}     
}


// Compute the function and its gradient
template<int dim>
void EigenvalueProblem<dim>::min_func_n_grad (const gsl_vector *x, double *f, gsl_vector *grad_func)
{
//try passing object by reference 
	*f  = min_func(x);
	min_func_grad(x , grad_func); 
}


//===============================================================================
//
//==============================================================================
//
//===============================================================================

template <int dim>
void EigenvalueProblem<dim>::assemble_system ()
{
   QGauss<dim>   quadrature_formula(8);

   FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values |
									update_gradients |
                           update_quadrature_points | 
									update_JxW_values);

	FEValues<dim> fe_values_poisson (fe_poisson, quadrature_formula,
									update_values |
									update_quadrature_points);

   const unsigned int dofs_per_cell = fe.dofs_per_cell;
   const unsigned int n_q_points    = quadrature_formula.size();

   FullMatrix<double> cell_stiffness_matrix (dofs_per_cell, dofs_per_cell);
   FullMatrix<double> cell_mass_matrix (dofs_per_cell, dofs_per_cell);
	
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);


	std::vector<double> potential_values (n_q_points);
	std::vector<double> cell_values_V_eff(n_q_points);	
	std::vector<double> cell_values_v_xc(n_q_points);

	
	typename DoFHandler<dim>:: active_cell_iterator
	cell_poisson = dof_handler_poisson.begin_active(),
	endc_poisson = dof_handler_poisson.end();

	typename DoFHandler<dim>::active_cell_iterator
 	cell = dof_handler.begin_active (),
	endc = dof_handler.end ();
//	pcout<<"v_xc size "<<v_xc.size() <<std::endl;
	if (XC_model_type == 0 )
	{
		if(scf_iteration==0)
		{
			for (; cell!=endc; ++cell)
	   	{
				if(cell->subdomain_id() == this_mpi_process)
				{			
	   	   	fe_values.reinit (cell);
	
					cell_stiffness_matrix = 0;
      			cell_mass_matrix = 0;
       			for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
					{	
						const Point<dim> pt = fe_values.quadrature_point(q_point);
	
         			for (unsigned int i=0; i<dofs_per_cell; ++i)
						{
            			for (unsigned int j=0; j<dofs_per_cell; ++j)
              			{
               			cell_stiffness_matrix (i, j)
						  		+= (0.5 * fe_values.shape_grad (i, q_point) *
   	               	 	 fe_values.shape_grad (j, q_point)
									 +
         	      	   	(v_external(fe_values.quadrature_point(q_point))/* + 
									optimum_v_kant.v_kantorovich(pt)*/) *
               	   		fe_values.shape_value (i, q_point) *
                 		   	fe_values.shape_value (j, q_point)
                  		  	) * fe_values.JxW (q_point);

		               	cell_mass_matrix (i, j)
   		   	          	+= (fe_values.shape_value (i, q_point) *
      		              	fe_values.shape_value (j, q_point)) *
									fe_values.JxW (q_point);
            	 		}
						}
					}	
					cell->get_dof_indices (local_dof_indices);
      		  	constraints.
					distribute_local_to_global (cell_stiffness_matrix,
   	      	                            local_dof_indices,
      	      	                         stiffness_matrix);
	        		constraints.
					distribute_local_to_global (cell_mass_matrix,
      	      	                         local_dof_indices,
         	      	                      mass_matrix);
//			 		++cell_poisson;
				}
			}
		}
		else
		{
 
			Potential_function<dim> optimum_v_kant(&optimum_NM_vars,
														no_electrons,
														nuclei_locations);
			const MappingQ1<dim> map_for_rho;
			std::vector<Point<dim>>  support_points_v(dof_handler.n_dofs());
//			std::map<types::global_dof_index, Point<dim>>  support_points_rho;
			DoFTools::map_dofs_to_support_points(map_for_rho,
															dof_handler, 
															support_points_v);
			solution_V_eff.reinit(dof_handler.n_dofs());	
			
			for(unsigned int i = 0; i< dof_handler.n_dofs(); ++i)
			{
				Point<3> origin = {0,0,0};
				//pcout<<"v_kanto at origin = "<< optimum_v_kant.v_kantorovich(origin)<<std::endl;
				const Point<3> sp = support_points_v[i];
				if(sp[0]==origin[0] & sp[1]==origin[0] & sp[2]==origin[2])
				{
					origin = {0.000001,0,0};
					solution_V_eff(i) = optimum_v_kant.v_kantorovich(origin);
				}
				else
				{
					solution_V_eff(i) = optimum_v_kant.v_kantorovich(sp);
				}

			}
//		constraints_rho.distribute(solution_V_eff);

			for (; cell!=endc; ++cell)
   		{
				if(cell->subdomain_id() ==this_mpi_process)
				{
					fe_values.reinit (cell);

					cell_stiffness_matrix = 0;
					cell_mass_matrix = 0;
					for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
					{	
						const Point<dim> pt = fe_values.quadrature_point(q_point);
		
						for (unsigned int i=0; i<dofs_per_cell; ++i)
						{
							for (unsigned int j=0; j<dofs_per_cell; ++j)
							{
								cell_stiffness_matrix (i, j)
								+= (0.5 * fe_values.shape_grad (i, q_point) *
									 fe_values.shape_grad (j, q_point)
									 +
									(v_external(fe_values.quadrature_point(q_point)) + 
									optimum_v_kant.v_kantorovich(pt)) *
									fe_values.shape_value (i, q_point) *
									fe_values.shape_value (j, q_point)
									) * fe_values.JxW (q_point);

								cell_mass_matrix (i, j)
								+= (fe_values.shape_value (i, q_point) *
									fe_values.shape_value (j, q_point)) *
									fe_values.JxW (q_point);
							}
						}
					}	
					cell->get_dof_indices (local_dof_indices);

					constraints.
					distribute_local_to_global (cell_stiffness_matrix,
														 local_dof_indices,
														 stiffness_matrix);
					constraints.
					distribute_local_to_global (cell_mass_matrix,
													  local_dof_indices,
													  mass_matrix);
	 			}	
			}
		}
	}	
	else
	{
		pcout<< "Error: The Exchange Correlation Model requested for is not implemented:	Select XC_model_type wihting range {0,1} "<< std::endl;
	}

	
    stiffness_matrix.compress (VectorOperation::add);
    mass_matrix.compress (VectorOperation::add);
}

template <int dim>
unsigned int EigenvalueProblem<dim>::solve ()
{
	double shift_constant = -20.0;  // value of shift in shift and invert transformation
	double max_iter_eigen_solver = 300;
	PETScWrappers::PreconditionerBase *preconditioner;
	preconditioner = new PETScWrappers::PreconditionBlockJacobi(mpi_communicator);
	SolverControl linear_solver_control (dof_handler.n_dofs(), 1e-15);
  
	PETScWrappers::SolverCG linear_solver(linear_solver_control,
			  										mpi_communicator);

	linear_solver.initialize(*preconditioner);

	SolverControl solver_control(max_iter_eigen_solver/*max no. of iterations*/,
											1e-15/*tol*/);

	SLEPcWrappers::SolverKrylovSchur *eigensolver;
	eigensolver = new SLEPcWrappers::SolverKrylovSchur(solver_control,mpi_communicator);
   SLEPcWrappers::TransformationShiftInvert::AdditionalData additional_data(shift_constant);
   SLEPcWrappers::TransformationShiftInvert 
						spectral_transformation(mpi_communicator, 
														additional_data);

   spectral_transformation.set_solver(linear_solver);
   eigensolver->set_transformation(spectral_transformation);

   eigensolver->set_which_eigenpairs (EPS_SMALLEST_REAL);

   eigensolver->set_problem_type (EPS_GHEP);

   eigensolver->solve (stiffness_matrix,
								mass_matrix,
                       	eigenvalues,
								eigenfunctions,
                       	eigenfunctions.size() );

// Eigenvectors should be normalized using the L2 norm.

	local_rho_current.reinit(dof_handler_poisson.n_dofs());
	Vector<double> local_temp_vec_rho;
	Vector<double> local_temp_eig_vector;
	double L2[eigenvalues.size()] = { };
	E_kin=0.0;	


	for (unsigned int i = 0; i<no_electrons / 2; ++i)
	{
     	Vector<double> local_temp_eig_vector(eigenfunctions[i]);
		constraints.distribute(local_temp_eig_vector);

		local_temp_vec_rho.reinit(dof_handler_poisson.n_dofs());
   	FETools::interpolate (dof_handler,
       	                  //constraints_hanging_eig_funs,
           	               local_temp_eig_vector,
              	            dof_handler_poisson,
                 	         constraints_rho,
                    	      local_temp_vec_rho);
		constraints_rho.distribute(local_temp_vec_rho);

		local_temp_vec_rho *= sqrt(2.0);

		PETScWrappers::MPI::Vector K_x(eigenfunctions[0]);
		stiffness_matrix.vmult(K_x,eigenfunctions[i]);
		E_kin += eigenfunctions[i] * K_x;


		QGauss<dim> quadrature_formula(8);	
		FEValues<dim> fe_values (fe_poisson, quadrature_formula,
										update_values |
										update_quadrature_points | 
										update_JxW_values);

	  	const unsigned int n_q_points = quadrature_formula.size();
		std::vector<double> cell_rho_current(n_q_points);

		typename DoFHandler<dim>::active_cell_iterator
		cell = dof_handler_poisson.begin_active (),
		endc = dof_handler_poisson.end ();
		for (; cell!=endc; ++cell)
		{
			 fe_values.reinit (cell);
			 fe_values.get_function_values(local_temp_vec_rho, cell_rho_current);
			 for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
			 {
				 L2[i] += cell_rho_current[q_point]*
							 cell_rho_current[q_point]*
							 fe_values.JxW(q_point);	 
			 }
		}

	   for (unsigned int j=0; j < dof_handler_poisson.n_dofs(); ++j)
      {
			eig_funcs_probability[i][j] = local_temp_vec_rho[j]*local_temp_vec_rho[j]; 
   		local_rho_current(j) += local_temp_vec_rho[j]*local_temp_vec_rho[j];
  	 	}
	
//		pcout<<"eige_vector :"<<i<<"integrates to \t"<<L2[i]<< std::endl;
	}
	delete preconditioner;
	delete eigensolver;	

	const double precision = 1e-5;
   PETScWrappers::MPI::Vector Kx(eigenfunctions[0]), Mx(eigenfunctions[0]), M2x(eigenfunctions[0]);
   for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
   {
   	mass_matrix.vmult(Mx, eigenfunctions[i]);
		for (unsigned int j = 0; j < eigenfunctions.size(); j++)
	   {
			double mass_norm = eigenfunctions[j] * Mx ; 
//			pcout	<<"mass_norm = "<< mass_norm
//					<<std::endl;
			Assert(std::abs(mass_norm-(i==j)) < precision,
                ExcMessage("Eigenvectors " + Utilities::int_to_string(i) +
   	                     " and " + Utilities::int_to_string(j) +
                           " are not orthonormal!"));
		}
		stiffness_matrix.vmult(Kx, eigenfunctions[i]);
      Kx.add(-1.0 * eigenvalues[i], Mx);
      Assert(Kx.l2_norm() < precision,
             ExcMessage(Utilities::to_string(Kx.l2_norm())));
    }

	return solver_control.last_step ();
}



template <int dim>
void EigenvalueProblem<dim>::E_xc_SCE()
{
	QGauss<dim> gauss_pts_per_dir(8); 
	FEValues<dim> cell_values(fe_poisson, gauss_pts_per_dir, 
									update_values | 
									update_quadrature_points | 
									update_JxW_values);	

	const unsigned int no_quad_points = gauss_pts_per_dir.size();

	typename DoFHandler<dim>::active_cell_iterator 
				cell = dof_handler_poisson.begin_active(),
				endc = dof_handler_poisson.end();

	std::vector<double> rho_next_cell (no_quad_points);
	std::vector<double> v_xc_cell (no_quad_points);

	Potential_function<dim> optimum_v_kant2(&optimum_NM_vars,
														no_electrons,
														nuclei_locations);

	E_xc = 0.0;
	for(;cell != endc; ++cell) 
	{
		cell_values.reinit(cell);
		
		cell_values.get_function_values (local_rho_current,rho_next_cell);

		for(unsigned int i=0; i < no_quad_points; ++i)
		{
			const Point<dim> pt = cell_values.quadrature_point(i);
			E_xc += rho_next_cell[i]*
					  optimum_v_kant2.v_kantorovich(pt) *
						cell_values.JxW(i);
		}
	}
	pcout<<"E_xc= \t"<<E_xc<<std::endl;	
	E_xchange_correlation.push_back(E_xc);
}


template <int dim>
void EigenvalueProblem<dim>::E_xc_n_Hartree()
{


	QGauss<dim> gauss_pts_per_dir(8); // this is nothing but the quadrature formula
	FEValues<dim> cell_values(fe, gauss_pts_per_dir, 
									update_values | 
									update_quadrature_points | 
									update_JxW_values);	

	FEValues<dim> fe_values_poisson(fe_poisson, gauss_pts_per_dir,
									update_values |
									update_quadrature_points);


	const unsigned int no_quad_points = gauss_pts_per_dir.size();

	typename DoFHandler<dim>::active_cell_iterator 
				cell = dof_handler.begin_active(),
				endc = dof_handler.end();
	typename DoFHandler<dim>::active_cell_iterator
				cell_poisson = dof_handler_poisson.begin_active(),
				endc_poisson = dof_handler_poisson.end();


	std::vector<double> rho_next_cell (no_quad_points);
	std::vector<double> solution_V_hartree_cell (no_quad_points);
	std::vector<double> v_xc_cell (no_quad_points);
	std::vector<double> epsilon_xc_cell (no_quad_points);


	E_H = 0.0;
	E_xc_n_H = 0.0;
	E_xc = 0.0;
	for(;cell != endc; ++cell) 
	{
		cell_values.reinit(cell);
		fe_values_poisson.reinit(cell_poisson);
		
		cell_values.get_function_values (epsilon_xc, epsilon_xc_cell);
		cell_values.get_function_values (v_xc, v_xc_cell);
		cell_values.get_function_values (local_rho_next,rho_next_cell);
		fe_values_poisson.get_function_values (solution_V_hartree,solution_V_hartree_cell);

		for(unsigned int i=0; i < no_quad_points; ++i)
		{
			E_xc_n_H += rho_next_cell[i]*
						(v_xc_cell[i] + 0.5*solution_V_hartree_cell[i] ) *
						cell_values.JxW(i);

			E_H += rho_next_cell[i] *
						0.5*solution_V_hartree_cell[i]*
						cell_values.JxW(i);	

			E_xc +=	rho_next_cell[i] *
						epsilon_xc_cell[i] * 
						cell_values.JxW(i);
		}
		++cell_poisson; 
	}
	pcout<<"Hartree energy= \t"<<E_H<<std::endl;	
	E_Hartree.push_back(E_H);
	E_xchange_correlation.push_back(E_xc);
	
}


template <int dim>
void EigenvalueProblem<dim>::energy()
{
	double sum_eigen_values = 0.0;
	double E_ground_state = 0.0;
	double f_i =2.0;
	
	for(int i=0; i < ceil(no_electrons/2.0); ++i)
	{
		if(i == floor(no_electrons/2.0))
			f_i =1.0;

		sum_eigen_values += f_i*eigenvalues[i];
	}
	if(XC_model_type==0)
	{
		if(scf_iteration==0)
		{		
			E_xchange_correlation.push_back(0.0);
			V_ee_SCE.push_back(0.0);
			E_ground_state = sum_eigen_values;
		}
		else
		{
			E_xc_SCE();
			E_ground_state = sum_eigen_values -
								E_xc - 
								V_ee_SCE[scf_iteration];
		}
		energy_values.push_back (E_ground_state);
	}
	else
	{	
		E_xc_n_Hartree();	
		E_ground_state = sum_eigen_values - 
							E_xc_n_H + 
							E_xc;
 
		energy_values.push_back (E_ground_state);
	}
	pcout<<"sum_eigen_values = "<<sum_eigen_values
				<<"\n E_ground_state =\t"<< E_ground_state
				<<"\nE_xc "<< E_xc
				<<"\n V_ee_SCE="<<V_ee_SCE[scf_iteration]
				<<std::endl;
	if(scf_iteration == 0)	
	{
		delta_E = 1.0;
	}
	else
	{
		delta_E = fabs(energy_values[scf_iteration-1] - energy_values[scf_iteration]);
	}
	pcout<<"Delta_e=\t"<<delta_E
				<<"E_ground_state=\t"<<E_ground_state<<std::endl;
}

//==================================================================================	

template <int dim>
void EigenvalueProblem<dim>::density_mixing()
{
	local_rho_next.reinit(dof_handler_poisson.n_dofs());
	double nu =0.0;
	if(scf_iteration < 5)
	{
		nu = 0.9; 
	}
	else
	{
		nu =0.9;
	}	

	if (scf_iteration==0)
	{
		for(unsigned int i=0; i < dof_handler_poisson.n_dofs(); ++i)
		{
			local_rho_next(i) = local_rho_current(i);
			local_rho_previous(i) = local_rho_next(i);
		}
	}
	else
	{	
		for(unsigned int i=0; i< dof_handler_poisson.n_dofs(); ++i)
		{
			local_rho_next(i)	= (1.0-nu) * local_rho_current(i) + nu * local_rho_previous(i);
			local_rho_previous(i) = local_rho_next(i);
		}
	}
}

//----------------End of Desnity mixing----------------------------------------
//=============================================================================




template <int dim>
void EigenvalueProblem<dim>::run ()
{
	if(XC_model_type == 0)
	{ 
		do
		{
			pcout<<"XC_model_type:SCE limit."<<std::endl
						<<"Begin scf_iteration loop #"<<scf_iteration<<": "<< std::endl;
			if (scf_iteration == 0)		
			{
				for (unsigned int cycle_apriori=0; cycle_apriori < refine_cycles; ++cycle_apriori)
				{  ref_cycle_no = cycle_apriori;
					pcout << "Refine_Cycle_apriori #" << cycle_apriori << ':' << std::endl;
      			if (cycle_apriori == 0)
					{
          			GridGenerator::hyper_cube (triangulation,-10,10);
          			triangulation.refine_global (3);
//						move_mesh();
						
        			}
      			else if (cycle_apriori >0 && cycle_apriori <2)
					{
						setup_system_V_hartree();
						refine_grid_apriori ();
					}
					else	
					{	
     					pcout << "   Number of active cells:       "
     		    			       << triangulation.n_active_cells()
     					          << std::endl;
					
						setup_system_V_hartree ();
      			
						pcout << "   Number of Poisson d.o.f: "
               				 << dof_handler_poisson.n_dofs()
			            	    << std::endl;
      				
//------- End of V_hartree refinement cycle for scf = 0 ----------
	 					setup_eigen_system ();

		    			assemble_system ();  
		    			const unsigned int n_iterations = solve ();
	   		 		pcout << "   Solver converged in " << n_iterations
         				   	 << " iterations." << std::endl;

                  if(cycle_apriori< refine_cycles-1)
                  {
                     refine_grid_apriori();
                  }
					}
				}
				density_mixing();
//				V_xc_SCE();
				energy();
//				output_results ();
				pcout << std::endl;
				for (unsigned int i=0; i<eigenvalues.size(); ++i)
      			pcout << "      Eigenvalue " << i
         			       << " : " << eigenvalues[i]
            	   		 << std::endl;

				pcout<<"############  SCF = 0 completed  ##############"<<std::endl;
			}
			else
			{
		//----------------------------
				setup_eigen_system ();
				V_xc_SCE();
			}

		pcout<<"End of SCF_iteration \t# "<<scf_iteration<<std::endl;

		if(delta_E <= E_tol)
		{
			break;
		}

		++scf_iteration;
		}
		while(delta_E >= E_tol && scf_iteration < 2);  //end of do loop
	}
} // end of run().


}// end of namespace Step-17

//=========================================================================
//=========================================================================
//=========================================================================


// main function:
int main (int argc, char **argv)
{
 
	int t0=time(NULL);
	std::cout<<"Trying"<<std::endl;
	try
   {
      using namespace dealii;
      using namespace Step36;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

		std::cout<<"Enter namespace Step36"<<std::endl;
      EigenvalueProblem<3> SCE ("step-36.prm");

      SCE.run ();
    }

    catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
		<< "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  std::cout << std::endl
            << "   Job done."
            << std::endl;

int tf = time(NULL);
std::cout<<"end time<<"<<tf<<std::endl;
printf ("Problem Completed! GSL used successfully! \n Time = %d seconds", tf-t0);
return 0;
}
