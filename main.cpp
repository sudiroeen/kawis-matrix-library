#include "kawis/kawis.hpp"

// using namespace kawis;

int main(){
	// Matrix ma(3,4);
	// cout << "ma: " << ma << endl;

	// Matrix m2 = {
	// 	{.1, .2, .3, .7},
	// 	{.4, .5, .6, .8},
	// 	{.9, .0, .1, .2}
	// };

	// std::unordered_map<std::string, Matrix> qrm2 = m2.qr();

	// cout << qrm2 << endl;

	// cout << "m2: " << m2 << endl;


	// concat({
	// 	{ma, ma, ma},
	// 	{m2, m2, m2}
	// });

	// Matrix m3 = m2;
	// cout << "m3: " << m3 << endl;

	// Matrix m4(ma.size());
	// cout << "m4: " << m4 << endl;
	// cout << "size m4: " << m4.size() << endl;

	// cout << "MatZeros(m4.T.size()): " << Matrix::Zeros(m4.t().size()) << endl;
	// cout << "MatZeros(5,4): " << Matrix::Zeros(5,4) << endl;

	// size_m sm5 = std::make_pair(2,7);
	// cout << "MatOnes" << sm5 << ": " << Matrix::Ones(sm5) << endl;

	// cout << "I3: " << Matrix::Identity(m2.row) << endl;

	// // operator* and friend operator* (float)
	// Matrix mn = float(2.) * Matrix::Ones(7,9) * 2.;
	// cout << "mn: " << mn << endl;

	// ma = Matrix::Ones(1,19);
	// cout << "ma: " << ma << endl;

	// // // operator=
	// Matrix mb = {
	// 	{1, 2, 3, 10},
	// 	{4, 5, 6, 11},
	// 	{7, 8, 9, 12},
	// 	{13, 14, 15, 16}
	// };

	// // operator* (Matrix)
	// ma = mb * Matrix::Ones(mb.col, 2);
	// cout << "ma: " << ma << endl;

	// // operator+ (Matrix)
	// ma = mb + Matrix::Ones(mb.size());
	// cout << "ma: " << ma << endl;

	// mb = {
	// 	{13, 14, 15, 16}
	// };
	// cout << "mb: " << mb << endl;

	// Matrix mv = {
	// 	{1, 2, 3, 10},
	// 	{4, 5, 6, 11},
	// 	{13, 14, 15, 16}
	// };

	// cout << "-mv: " << -mv << endl;
	// cout << "mv-m2: " << mv-m2 << endl;
	// cout << "mv/5.: " << mv/5. << endl;
	// cout << "mv^T: " << mv.t() << endl;

	// Matrix m7 = Matrix::Ones(mb.col, 1) * mb;
	// cout << "m7: " << m7 << endl;
	// cout << "m7D: " << m7.getDiag() << endl;

	// m7.setDiag({{0, 0, 0, 0}});
	// cout << "m7New: " << m7 << endl;

	// Matrix mvBlock = mv.block(0,0,1,3);
	// cout << "mvBlock: " << mvBlock << endl;

	// mv.setBlock(0,0,1,2, 999.*Matrix::Ones(2,3));
	// cout << "mvNewBlock: " << mv << endl;

	// mv.setBlock(0, 0, 0, 2, {{55, 55, 55}});
	// cout << "mvNewBlock: " << mv << endl;

	// cout << "size of mvBlock: " << mvBlock.size() << endl;
	
	// Matrix rescon = concat(
	// 	{
	// 		{10.*m2, Matrix::Zeros(3,1), Matrix::Identity(3)},
	// 		{Matrix::Identity(5), Matrix::Ones(5,3)        }
	// 	});

	// rescon.print("rescon");

	// cout << "det: " << rescon.det() << endl;

	// Matrix tath = {
	// 	{0, 0, 0, 0, 2, 0, 3},
	// 	{0, 0, 3, 0, 4, 0, 4},
	// 	{0, 0, 5, 0, 6, 0, 6},
	// 	{0, 0, 7, 0, 8, 0, 8},
	// 	{0, 0, 9, 0, 10, 0, 10},
	// 	{0, 0, 10, 0, 11, 0, 11},
	// 	{0, 0, 11, 0, 12, 0, 12}
	// };

	// std::unordered_map<std::string, Matrix> tathqr = tath.qr();
	// cout << tathqr << endl;

	// Matrix Q = tathqr["Q"];
	// Matrix R = tathqr["R"];

	// (Q * R).print("QR");

	// cout << "is same tath: " << (Q * R).isSame(tath) << endl;

	// for(int c=0; c<Q.col; c++)
	// 	cout << Q.getCol(c).norm() << endl;

	// cout << "is same: " << (Q * Q.t()).isSame(Matrix::Identity(Q.row)) << endl;

	// tath.deleteCol(0);

	// std::unordered_map<std::string, Matrix> LUrescon = rescon.LUdecomp();
	// cout << LUrescon << endl;

	// cout << "E: " << Emat << endl;
	// cout << "L: " << Lmat << endl;
	// cout << "U: " << Umat << endl;

	// P*A = L*U -> A = Pinv*L*U = Ptrans * L * U
	// cout << "P*L*U: " << Pmat.t() * Lmat * Umat << endl;

	// L = Einv and E = Linv then E*L = L*E = I
	// cout << "L*E: " << Lmat*Emat << "\nE*L: " << Emat*Lmat << endl;

	// std::unordered_map<std::string, Matrix> LUZero = Matrix::Zeros(2,3).LUdecomp();
	// cout << LUZero << endl;

	// tath.print("tath");
	// std::unordered_map<std::string, Matrix> LUtath = tath.LUdecomp();
	// cout << LUtath << endl;

	// Matrix tath_t = tath.t();
	// tath_t.print("tath_t");
	// std::unordered_map<std::string, Matrix> LUtath_t = tath_t.LUdecomp();
	// cout << LUtath_t << endl;

	// Matrix P = LUtath_t["P"];
	// Matrix L = LUtath_t["L"];
	// Matrix E = LUtath_t["E"];
	// Matrix U = LUtath_t["U"];

	// P*A = L*U
	// A = Pinv * L * U = Ptrans * L * U
	// Linv * P * A = U sehingga E * P * A = U
	// (P.t() * L * U).print("P.t * L * U");
	// (E * P * tath.t()).print("E * P * A");

	// tath.getDiag().print("tathD");
	// tath.getRow(2,-1).getDiag().print("tath_2_-1_D");

	// cout << rescon.det(1) << endl;
	// LUrescon["U"].getDiag().print("diagU");

	// cout << tath.t().rank() << endl;

	// cout << rescon.rank() << endl;

	// std::unordered_map<std::string, Matrix> res_rref = rescon.rref();
	// res_rref["R"].print("RREF");

	// Matrix invA = res_rref["E"];
	// (invA * rescon).print("is I");
	// (rescon * invA).print("is I");

	// std::unordered_map<std::string, Matrix> tath_rref = tath.t().rref();
	// tath_rref["R"].print("RREF TATH");

	// Matrix m3 = {
	// 	{1, 2, 3, 5, 7},
	// 	{3, 6, 9, 13, 17}, //{2, 4, 6, 12, 7},
	// 	{3, 6, 9, 13, 18}
	// 	// {0, 0, 0, 1, 1} //{3, 6, 9, 13, 17}
	// };

	// std::unordered_map<std::string, Matrix> m3_lu = m3.LUdecomp();
	// cout << m3_lu << endl;

	// (m3_lu["L"] * m3_lu["U"] * m3_lu["PR"]).print("LU");

	// (m3_lu["U"] * m3_lu["PR"]).print("U*PR");

	// std::unordered_map<std::string, Matrix> m3_rref = m3.rref();

	// m3_rref["rref"].print("m3_rref");

	// Matrix ns_m3 = m3.null_space();

	// ns_m3.print("Nullspace m3");

	// Matrix Amat = {
	// 	{1, 2, 1},
	// 	{3, 8, 3}, //{3, 8, 1},
	// 	{4, 8, 4}//{0, 1, 1}
	// };

	// Matrix bmat = {
	// 	{3, 2, 7},
	// 	{11, 12, 21},
	// 	{12, 2, 28}
	// };

	// Amat.print("Amat");
	// bmat.print("bmat");

	// std::unordered_map<std::string, Matrix> xsolve = solve(Amat, bmat);
	// cout << "xsolve: \n" << xsolve << endl;

	// float a = std::numeric_limits<float>::quiet_NaN();
	// cout << "a: " << a << endl;

	// Matrix::MatrixNaN(5,3).print("NaN matrix");

	// Matrix A = {
	// 	{1,2,3},
	// 	{4,5,6},
	// 	{7,8,10}
	// };
	// std::unordered_map<std::string, Matrix> res_qr_iter = A.qr_iter(EPSILON, -1, false);

	// cout << A.eig() << endl;

	Matrix Amat = {
		{20, 1, 1, 0, 0, 0},
		{5, -2, 0, -1, 0, 0},
		{4, -1, 0, 0, -1, 0},
		{0, 3, 0, 0, 0, -1},
		{0, 0, -1, 0, 1, 0},
		{0, 0, 0, -1, 0, 1}
	};

	Matrix bmat = {
		{40},
		{-16},
		{-4},
		{4},
		{0},
		{0}
	};

	// concat({{Amat, bmat}}).LUdecomp();

	// concat({{Amat, bmat}}).rref();

	// cout << "Ludecomp: " << Amat.LUdecomp() << endl;

	cout << "solveAb: " << solve(Amat, bmat) << endl;

	// Matrix A1 = {
	// 	{-1, 1, -1},
	// 	{16, 0, 12},
	// 	{12, 0, 20}
	// };

	// Matrix b1 = {
	// 	{0},
	// 	{36},
	// 	{24}
	// };

	// cout << "solve I: " << solve(A1, b1*11.) << endl;

	return 0;
}