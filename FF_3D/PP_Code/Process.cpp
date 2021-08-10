//#define PRINT_DEBUGGING_OUTPUT
#define USE_SCALAR_IMPLEMENTATION

#define COMPUTE_V_AS_MATRIX
//#define COMPUTE_V_AS_QUATERNION
#define COMPUTE_U_AS_MATRIX
//#define COMPUTE_U_AS_QUATERNION

#define min_value 1e-20

#include "Singular-Value-Decomposition-master\Singular_Value_Decomposition_Preamble.hpp"
#include "Process.h"
#include "OMPHelper.h"

Process::Process(string filename, string init_filename)
{
	fileManager.readFile(filename, source_mesh);
	fileManager.readFile(init_filename, init_mesh);
	
	lb = -1;
	ub = -1;
	iter_max = 1000;
	tol_err = 1e-4;
	tol_energy = 1e-10;
	total_iter = 0;

	vector<double> & cell_vol = source_mesh.volume;
	cell_vol.resize(source_mesh.n_cells());
	OpenVolumeMesh::CellIter c_it = source_mesh.cells_begin();
	OpenVolumeMesh::CellVertexIter cv_it = source_mesh.cv_iter(*c_it);
	Eigen::Matrix4d volume_matrix;
	OpenVolumeMesh::Geometry::Vec3d P; int count = 0;
	for (c_it; c_it != source_mesh.cells_end(); ++c_it)
	{
		cv_it = source_mesh.cv_iter(*c_it);
		count = 0;
		for (cv_it; cv_it.valid(); ++cv_it)
		{
			P = source_mesh.vertex(*cv_it);
			volume_matrix(1, count) = P[0];
			volume_matrix(2, count) = P[1];
			volume_matrix(3, count) = P[2];
			volume_matrix(0, count) = 1.0;
			++count;
		}

		cell_vol[c_it->idx()] = volume_matrix.determinant() / 6.0;
	}
}

Process::~Process()
{
	if (pardiso != NULL)
	{
		delete pardiso;
		pardiso = NULL;
	}
}

bool Process::writeobj(VolumeMesh& out_mesh, const string out_filename)
{
	for (auto it1 = out_mesh.vertices_begin(); it1 != out_mesh.vertices_end(); ++it1)
	{
		int id = it1->idx();
		OpenVolumeMesh::Geometry::Vec3d pos(cur_position(id), cur_position(id + V_N), cur_position(id + 2 * V_N));
		out_mesh.set_vertex(*it1, pos);
	}
	return fileManager.writeFile(out_filename, out_mesh);
}

bool Process::writeobj(const string out_filename)
{
	for (auto it1 = init_mesh.vertices_begin(); it1 != init_mesh.vertices_end(); ++it1)
	{
		int id = it1->idx();
		OpenVolumeMesh::Geometry::Vec3d pos(cur_position(id), cur_position(id + V_N), cur_position(id + 2 * V_N));
		init_mesh.set_vertex(*it1, pos);
	}
	return fileManager.writeFile(out_filename, init_mesh);
}

void Process::init()
{
	V_N = source_mesh.n_vertices();
	F_N = source_mesh.n_cells();

	F.resize(F_N, 4);
	OpenVolumeMesh::CellIter  f_it = source_mesh.cells_begin();
	OpenVolumeMesh::CellIter f_end = source_mesh.cells_end();
	OpenVolumeMesh::CellVertexIter fv_it = source_mesh.cv_iter(*f_it);

	source_position.resize(V_N * 3);
	OpenVolumeMesh::VertexIter v_it = source_mesh.vertices_begin();
	OpenVolumeMesh::VertexIter v_end = source_mesh.vertices_end();
	VolumeMesh::PointT pos;

	global_position.resize(V_N * 3);
	OpenVolumeMesh::VertexIter iv_it = init_mesh.vertices_begin();
	VolumeMesh::PointT ipos;

	Tx.resize(F_N * 3 * 3);
	pTx.resize(F_N * 3 * 3);
	distortions.resize(F_N);
	flips.resize(F_N);
	K.resize(F_N);

	int i = 0, k = 0;
	int T_rows = F_N * 3 * 3;
	int T_cols = V_N * 3;
	T.resize(T_cols, T_rows);
	vector<Tri> triplist_T;
	triplist_T.reserve(T_rows * 4);
	vector<Tri> triplist_A;
	triplist_A.reserve(2 * T_cols);

	double x12, x13, x14, y12, y13, y14, z12, z13, z14, a11, a12, a13, a21, a22, a23, a31, a32, a33, vol;

	MatrixXd B = MatrixXd::Identity(4, 4);
	B = B.array() - (1.0 / 4);

	T.resize(T_cols, T_rows);
	MatrixXd currV(4, 3);
	MatrixXd currT(3, 4);
	int curr_row = 0;
	n_constraints = 0;

	OMP_PARALLEL
	{
		OMP_FOR
		for (; f_it != f_end; ++f_it, ++i)
		{
			fv_it = source_mesh.cv_iter(*f_it);
			for (int j = 0; j < 4; ++j, ++fv_it)
			{
				F(i, j) = fv_it.cur_handle();
			}
		}

		OMP_FOR
		for (; v_it != v_end; ++v_it, ++k, ++iv_it)
		{
			pos = source_mesh.vertex(*v_it);
			ipos = init_mesh.vertex(*iv_it);
			for (int j = 0; j < 3; ++j)
			{
				source_position[k + V_N*j] = pos[j];
				global_position[k + V_N*j] = ipos[j];
			}
		}

		OMP_FOR
		for (auto it1 = init_mesh.vertices_begin(); it1 != init_mesh.vertices_end(); ++it1)
		{
			if (init_mesh.is_boundary(*it1))
			{
				triplist_A.push_back(Tri(3 * n_constraints, it1->idx(), 1));
				triplist_A.push_back(Tri(3 * n_constraints + 1, it1->idx() + V_N, 1));
				triplist_A.push_back(Tri(3 * n_constraints + 2, it1->idx() + 2 * V_N, 1));
				triplist_A.push_back(Tri(3 * n_constraints, T_cols + 3 * n_constraints, 0.0));
				triplist_A.push_back(Tri(3 * n_constraints + 1, T_cols + 3 * n_constraints + 1, 0.0));
				triplist_A.push_back(Tri(3 * n_constraints + 2, T_cols + 3 * n_constraints + 2, 0.0));
				n_constraints++;
			}
		}

		OMP_FOR
		for (int ii = 0; ii < F_N; ii++)
		{
			vol = 6 * source_mesh.volume[ii];
			x12 = source_position[F(ii, 1)] - source_position[F(ii, 0)];
			x13 = source_position[F(ii, 2)] - source_position[F(ii, 0)];
			x14 = source_position[F(ii, 3)] - source_position[F(ii, 0)];
			y12 = source_position[F(ii, 1) + V_N] - source_position[F(ii, 0) + V_N];
			y13 = source_position[F(ii, 2) + V_N] - source_position[F(ii, 0) + V_N];
			y14 = source_position[F(ii, 3) + V_N] - source_position[F(ii, 0) + V_N];
			z12 = source_position[F(ii, 1) + 2 * V_N] - source_position[F(ii, 0) + 2 * V_N];
			z13 = source_position[F(ii, 2) + 2 * V_N] - source_position[F(ii, 0) + 2 * V_N];
			z14 = source_position[F(ii, 3) + 2 * V_N] - source_position[F(ii, 0) + 2 * V_N];

			a11 = y13*z14 - y14*z13;	a12 = y14*z12 - y12*z14;	a13 = y12*z13 - y13*z12;
			a21 = x14*z13 - x13*z14;	a22 = x12*z14 - x14*z12;	a23 = x13*z12 - x12*z13;
			a31 = x13*y14 - x14*y13;	a32 = x14*y12 - x12*y14;	a33 = x12*y13 - x13*y12;

			triplist_T.push_back(Tri(F(ii, 0), curr_row, -(a11 + a12 + a13) / vol));
			triplist_T.push_back(Tri(F(ii, 1), curr_row, a11 / vol));
			triplist_T.push_back(Tri(F(ii, 2), curr_row, a12 / vol));
			triplist_T.push_back(Tri(F(ii, 3), curr_row, a13 / vol)); curr_row = curr_row + 1;
			triplist_T.push_back(Tri(F(ii, 0) + V_N, curr_row, -(a11 + a12 + a13) / vol));
			triplist_T.push_back(Tri(F(ii, 1) + V_N, curr_row, a11 / vol));
			triplist_T.push_back(Tri(F(ii, 2) + V_N, curr_row, a12 / vol));
			triplist_T.push_back(Tri(F(ii, 3) + V_N, curr_row, a13 / vol)); curr_row = curr_row + 1;
			triplist_T.push_back(Tri(F(ii, 0) + 2 * V_N, curr_row, -(a11 + a12 + a13) / vol));
			triplist_T.push_back(Tri(F(ii, 1) + 2 * V_N, curr_row, a11 / vol));
			triplist_T.push_back(Tri(F(ii, 2) + 2 * V_N, curr_row, a12 / vol));
			triplist_T.push_back(Tri(F(ii, 3) + 2 * V_N, curr_row, a13 / vol)); curr_row = curr_row + 1;

			triplist_T.push_back(Tri(F(ii, 0), curr_row, -(a21 + a22 + a23) / vol));
			triplist_T.push_back(Tri(F(ii, 1), curr_row, a21 / vol));
			triplist_T.push_back(Tri(F(ii, 2), curr_row, a22 / vol));
			triplist_T.push_back(Tri(F(ii, 3), curr_row, a23 / vol)); curr_row = curr_row + 1;
			triplist_T.push_back(Tri(F(ii, 0) + V_N, curr_row, -(a21 + a22 + a23) / vol));
			triplist_T.push_back(Tri(F(ii, 1) + V_N, curr_row, a21 / vol));
			triplist_T.push_back(Tri(F(ii, 2) + V_N, curr_row, a22 / vol));
			triplist_T.push_back(Tri(F(ii, 3) + V_N, curr_row, a23 / vol)); curr_row = curr_row + 1;
			triplist_T.push_back(Tri(F(ii, 0) + 2 * V_N, curr_row, -(a21 + a22 + a23) / vol));
			triplist_T.push_back(Tri(F(ii, 1) + 2 * V_N, curr_row, a21 / vol));
			triplist_T.push_back(Tri(F(ii, 2) + 2 * V_N, curr_row, a22 / vol));
			triplist_T.push_back(Tri(F(ii, 3) + 2 * V_N, curr_row, a23 / vol)); curr_row = curr_row + 1;

			triplist_T.push_back(Tri(F(ii, 0), curr_row, -(a31 + a32 + a33) / vol));
			triplist_T.push_back(Tri(F(ii, 1), curr_row, a31 / vol));
			triplist_T.push_back(Tri(F(ii, 2), curr_row, a32 / vol));
			triplist_T.push_back(Tri(F(ii, 3), curr_row, a33 / vol)); curr_row = curr_row + 1;
			triplist_T.push_back(Tri(F(ii, 0) + V_N, curr_row, -(a31 + a32 + a33) / vol));
			triplist_T.push_back(Tri(F(ii, 1) + V_N, curr_row, a31 / vol));
			triplist_T.push_back(Tri(F(ii, 2) + V_N, curr_row, a32 / vol));
			triplist_T.push_back(Tri(F(ii, 3) + V_N, curr_row, a33 / vol)); curr_row = curr_row + 1;
			triplist_T.push_back(Tri(F(ii, 0) + 2 * V_N, curr_row, -(a31 + a32 + a33) / vol));
			triplist_T.push_back(Tri(F(ii, 1) + 2 * V_N, curr_row, a31 / vol));
			triplist_T.push_back(Tri(F(ii, 2) + 2 * V_N, curr_row, a32 / vol));
			triplist_T.push_back(Tri(F(ii, 3) + 2 * V_N, curr_row, a33 / vol)); curr_row = curr_row + 1;
		}
	}

	SparseMatrix<double> A_extension(3 * n_constraints, T_cols + 3 * n_constraints);
	A_extension.setFromTriplets(triplist_A.begin(), triplist_A.end());
	b.resize(3 * n_constraints);
	for (int i = 0; i < n_constraints; ++i)
	{
		b(3 * i) = global_position(triplist_A[6 * i].col());
		b(3 * i + 1) = global_position(triplist_A[6 * i + 1].col());
		b(3 * i + 2) = global_position(triplist_A[6 * i + 2].col());
	}

	T.setFromTriplets(triplist_T.begin(), triplist_T.end());
	T = T.transpose();
	SparseMatrix<double> Coef;
	SparseMatrix<double> T_extension(T_cols + 3 * n_constraints, T_rows);
	T_extension.setFromTriplets(triplist_T.begin(), triplist_T.end());
	Coef = T_extension * T_extension.transpose();
	Coef.rightCols(3 * n_constraints) = A_extension.transpose();
	Coef = Coef.triangularView<Upper>();
	Coef = Coef.transpose();

	pardiso = new PardisoSolver();

	int * ia = Coef.outerIndexPtr();
	vector<int> pardiso_ia(ia, ia + 3 * V_N + 1 + 3 * n_constraints);
	pardiso->ia = pardiso_ia;

	int * ja = Coef.innerIndexPtr();
	int nonzero = Coef.nonZeros();
	vector<int> pardiso_ja(ja, ja + nonzero);
	pardiso->ja = pardiso_ja;

	double * a = Coef.valuePtr();
	vector<double> pardiso_a(a, a + nonzero);
	pardiso->a = pardiso_a;

	pardiso->nnz = pardiso_ja.size();
	pardiso->num = V_N * 3 + 3 * n_constraints;

	pardiso->pardiso_init();

	pardiso->factorize();
}

void Process::local()
{
	pTx = Tx;
	int block_size = 3 * 3;
	int num_blocks = pTx.size() / block_size;
	double K2plus1;
	double K2plus2;
	double TwoK2plus1;
	JacobiSVD<Matrix3d> svdA(3, 3, (ComputeFullU | ComputeFullV));
	Vector3d s, snew;
	double dist, t;
	bool flipped;
	Matrix3d U, V, U1, V1;
	Map<Matrix3d> currA(pTx.data());
	bool wasProjected;

	int c = 0;
	OMP_PARALLEL
	{
		OMP_FOR
		for (int ii = 0; ii < num_blocks; ii++)
		{
			K2plus1 = K[ii] * K[ii] + 1;
			K2plus2 = K[ii] * K[ii] + 2;
			TwoK2plus1 = 2 * K[ii] * K[ii] + 1;
			new (&currA) Map<Matrix3d>(pTx.data() + ii*block_size);

			if (abs(currA.determinant()) < min_value)
			{
				flipped = true;
				flips(ii) = flipped;
				distortions[ii] = std::numeric_limits<double>::infinity();
				currA.setIdentity();
				currA(0, 0) = K[ii] * min_value;
				currA(1, 1) = min_value;
				currA(2, 2) = min_value;
				continue;
			}

			flipped = (currA.determinant() < 0);
			flips(ii) = flipped;

#include "Singular-Value-Decomposition-master\Singular_Value_Decomposition_Kernel_Declarations.hpp"
			ENABLE_SCALAR_IMPLEMENTATION(Sa11.f = currA(0, 0);)
				ENABLE_SCALAR_IMPLEMENTATION(Sa21.f = currA(1, 0);)
				ENABLE_SCALAR_IMPLEMENTATION(Sa31.f = currA(2, 0);)
				ENABLE_SCALAR_IMPLEMENTATION(Sa12.f = currA(0, 1);)
				ENABLE_SCALAR_IMPLEMENTATION(Sa22.f = currA(1, 1);)
				ENABLE_SCALAR_IMPLEMENTATION(Sa32.f = currA(2, 1);)
				ENABLE_SCALAR_IMPLEMENTATION(Sa13.f = currA(0, 2);)
				ENABLE_SCALAR_IMPLEMENTATION(Sa23.f = currA(1, 2);)
				ENABLE_SCALAR_IMPLEMENTATION(Sa33.f = currA(2, 2);)

#include "Singular-Value-Decomposition-master\Singular_Value_Decomposition_Main_Kernel_Body.hpp"
				ENABLE_SCALAR_IMPLEMENTATION(U(0, 0) = Su11.f;)
				ENABLE_SCALAR_IMPLEMENTATION(U(1, 0) = Su21.f;)
				ENABLE_SCALAR_IMPLEMENTATION(U(2, 0) = Su31.f;)
				ENABLE_SCALAR_IMPLEMENTATION(U(0, 1) = Su12.f;)
				ENABLE_SCALAR_IMPLEMENTATION(U(1, 1) = Su22.f;)
				ENABLE_SCALAR_IMPLEMENTATION(U(2, 1) = Su32.f;)
				ENABLE_SCALAR_IMPLEMENTATION(U(0, 2) = Su13.f;)
				ENABLE_SCALAR_IMPLEMENTATION(U(1, 2) = Su23.f;)
				ENABLE_SCALAR_IMPLEMENTATION(U(2, 2) = Su33.f;)

				ENABLE_SCALAR_IMPLEMENTATION(V(0, 0) = Sv11.f;)
				ENABLE_SCALAR_IMPLEMENTATION(V(1, 0) = Sv21.f;)
				ENABLE_SCALAR_IMPLEMENTATION(V(2, 0) = Sv31.f;)
				ENABLE_SCALAR_IMPLEMENTATION(V(0, 1) = Sv12.f;)
				ENABLE_SCALAR_IMPLEMENTATION(V(1, 1) = Sv22.f;)
				ENABLE_SCALAR_IMPLEMENTATION(V(2, 1) = Sv32.f;)
				ENABLE_SCALAR_IMPLEMENTATION(V(0, 2) = Sv13.f;)
				ENABLE_SCALAR_IMPLEMENTATION(V(1, 2) = Sv23.f;)
				ENABLE_SCALAR_IMPLEMENTATION(V(2, 2) = Sv33.f;)

				ENABLE_SCALAR_IMPLEMENTATION(s[0] = Sa11.f;)
				ENABLE_SCALAR_IMPLEMENTATION(s[1] = Sa22.f;)
				ENABLE_SCALAR_IMPLEMENTATION(s[2] = Sa33.f;)

				dist = s(0) / s(2);
			if (flipped)
			{
				dist = -dist;
			}

			wasProjected = false;
			snew = s;
			if ((K[ii] > 0) & (flipped | (dist > K[ii])))
			{
				t = (K[ii] * s(0) + s(2)) / K2plus1;
				snew << t*K[ii], s(1), t;
				if (snew(1)<snew(2))
				{
					t = (K[ii] * s(0) + s(1) + s(2)) / K2plus2;
					snew << t*K[ii], t, t;
				}
				else if (snew(1)>snew(0))
				{
					t = (K[ii] * (s(0) + s(1)) + s(2)) / TwoK2plus1;
					snew << t*K[ii], t*K[ii], t;
				}
				wasProjected = true;
			}
			if ((lb > 0) & (snew(2) < lb))
			{
				snew = snew.array().max(lb);
				wasProjected = true;
			}
			if ((ub > 0) & (snew(0) > ub))
			{
				snew = snew.array().min(ub);
				wasProjected = true;
			}
			if (wasProjected)
			{
				currA = U*snew.asDiagonal()*V.transpose();
			}

			distortions(ii) = dist;
		}
	}
}

double Process::Energy()
{
	VectorXd energy1 = Tx - pTx;
	return energy1.norm();
}

void Process::pck()
{
	init();
	int iter = 0;
	cur_position = global_position;		//Q_0
	Tx = T*cur_position;
	Distortion();

	long time_beg, time_end;
	time_beg = clock();

	while (flips.sum() && iter< 10)
	{
		std::cout << "--------------------------- " << to_string(iter) << " times projection" << " ---------------------------" << std::endl;
		std::cout << "max distortions:   " << distortions.maxCoeff() << std::endl;
		std::cout << "flip number:       " << flips.sum() << std::endl;

		// update conformal distortion bound K (Section 3.4)
		for (int i = 0; i < F_N; ++i)
		{
			K(i) = 4 * pow(2, iter);
		}

		// project the mapping into the bounded distortion space (Section 3.3)
		solve();
		++iter;
	}

	time_end = clock();
	double time_consumption = (time_end - time_beg) / 1000.0;

	std::cout << "=========================== " << "optimal result data" << " ===========================" << std::endl;
	std::cout << "max distortions:       " << distortions.maxCoeff() << std::endl;
	double avg_distortion = 0;
	for (int i = 0; i < F_N; ++i)
	{
		avg_distortion += distortions(i);
	}
	avg_distortion = avg_distortion / F_N;
	std::cout << "average distortions:   " << avg_distortion << std::endl;
	std::cout << "flip number:           " << flips.sum() << std::endl;
	std::cout << "time_consumption:      " << time_consumption << std::endl;
	std::cout << "number of rounds:      " << iter << std::endl;
	std::cout << "number of solving equ: " << total_iter << std::endl;
	writeobj(init_mesh, "result.ovm");
}

void Process::solve()
{
	Anderson accelerator;
	int Anderson_m = 5;
	accelerator.init(Anderson_m, 3, cur_position);
	pre_energy = std::numeric_limits<double>::max();
	double is_convergence;

	for (int iter = 0; iter < iter_max; ++iter)
	{
		Tx = T*cur_position;
		local();

		if (Energy()>pre_energy)
		{
			cur_position = global_position;
			Tx = T*cur_position;
			local();

			// guarante to decrease the target energy at each iteration
			accelerator.replace(cur_position);
		}

		cur_energy = Energy();
		std::cout << iter << "    " << cur_energy << "    " << flips.sum() << std::endl;

		is_convergence = abs(pre_energy - cur_energy) / pre_energy;
		if (is_convergence < tol_err || (0 == flips.sum()) || cur_energy<tol_energy)
		{
			break;
		}
		pre_energy = cur_energy;

		// Anderson acceleration
		if ((iter_max - 1) != iter)
		{
			global();				
			cur_position = accelerator.compute(global_position);
		}
	}
}

void Process::global()
{
	VectorXd right(3 * V_N + 3 * n_constraints);
	right.head(3 * V_N) = T.transpose()*pTx;
	right.bottomRows(3 * n_constraints) = b;

	vector<double> pardiso_u(right.data(), right.data() + 3 * V_N + 3 * n_constraints);
	pardiso->rhs = pardiso_u;
	pardiso->pardiso_solver();

	for (size_t i = 0; i < 3 * V_N; i++)
	{
		global_position(i) = (pardiso->result)[i];
	}

	total_iter++;
}

void Process::Distortion()
{
	int block_size = 3 * 3;
	int num_blocks = Tx.size() / block_size;
	Vector3d s;
	double dist, t;
	bool flipped;
	Map<Matrix3d> currA(Tx.data());

	for (int ii = 0; ii < num_blocks; ii++)
	{
		new (&currA) Map<Matrix3d>(Tx.data() + ii*block_size);

		if (abs(currA.determinant()) < min_value)
		{
			flipped = true;
			flips(ii) = flipped;
			distortions[ii] = std::numeric_limits<double>::infinity();
			continue;
		}

		flipped = (currA.determinant() < 0);
		flips(ii) = flipped;

#include "Singular-Value-Decomposition-master\Singular_Value_Decomposition_Kernel_Declarations.hpp"
		ENABLE_SCALAR_IMPLEMENTATION(Sa11.f = currA(0, 0);)
			ENABLE_SCALAR_IMPLEMENTATION(Sa21.f = currA(1, 0);)
			ENABLE_SCALAR_IMPLEMENTATION(Sa31.f = currA(2, 0);)
			ENABLE_SCALAR_IMPLEMENTATION(Sa12.f = currA(0, 1);)
			ENABLE_SCALAR_IMPLEMENTATION(Sa22.f = currA(1, 1);)
			ENABLE_SCALAR_IMPLEMENTATION(Sa32.f = currA(2, 1);)
			ENABLE_SCALAR_IMPLEMENTATION(Sa13.f = currA(0, 2);)
			ENABLE_SCALAR_IMPLEMENTATION(Sa23.f = currA(1, 2);)
			ENABLE_SCALAR_IMPLEMENTATION(Sa33.f = currA(2, 2);)

#include "Singular-Value-Decomposition-master\Singular_Value_Decomposition_Main_Kernel_Body.hpp"

			ENABLE_SCALAR_IMPLEMENTATION(s[0] = Sa11.f;)
			ENABLE_SCALAR_IMPLEMENTATION(s[1] = Sa22.f;)
			ENABLE_SCALAR_IMPLEMENTATION(s[2] = Sa33.f;)

			dist = s(0) / s(2);
			distortions[ii] = dist;
	}
}
