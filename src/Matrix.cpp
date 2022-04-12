/*
	Copyright:
				Sudiro
						[at] SudiroEEN@gmail.com

*/
#include "../include/Matrix.hpp"

/////////////////////////  CONSTRUCTOR
Matrix::Matrix(){
	row = 0;
	col = 0;
}

Matrix::Matrix(int row_, int col_)
	   : row(row_), col(col_)
{
	entry = new double*[row];
	for(int r=0; r<row; r++){
		entry[r] = new double[col];
		for(int c=0; c<col; c++){
			entry[r][c] = 0.;
		}
	}
}

Matrix::Matrix(size_m sm)
	   : row(sm.first), col(sm.second)
{
	entry = new double*[row];
	for(int r=0; r<row; r++){
		entry[r] = new double[col];
		for(int c=0; c<col; c++)
			entry[r][c] = double(0.);
	}
}

Matrix::Matrix(const initializer_list< initializer_list <double> >& iif){
	row = (int)iif.size();
	col = (int)iif.begin()->size();

	entry = new double*[row];
	
	for(int r=0; r<row; r++){
		entry[r] = new double[col];
		for(int c=0; c<col; c++){
			entry[r][c] = ((iif.begin()+r)->begin())[c];
		}
	}
}

Matrix::Matrix(const Matrix& M){
	row = M.row;
	col = M.col;
	
	entry = new double*[row];
	for(int r=0; r<row; r++){
		entry[r] = new double[col];
		for(int c=0; c<col; c++)
			entry[r][c] = M.entry[r][c];
	}
}

//////////////////////// MEMBER FUNCTION BIASA
Matrix Matrix::t(){
	Matrix trp(col, row);
	for(int rt=0; rt<col; rt++){
		for(int ct=0; ct<row; ct++)
			trp.set(rt,ct) = (*this)(ct,rt);
	}
	return trp;
}

Matrix Matrix::getDiag(){
	int roc = n_diag();
	Matrix mdiag(roc, 1);
	for(int d=0; d<roc; d++)
		mdiag.set(d) = (*this)(d,d);
	return mdiag;
}

void Matrix::setDiag(const Matrix& rowmat){
	for(int d=0; d<rowmat.row; d++)
		entry[d][d] = rowmat.entry[d][0];
}

Matrix Matrix::block(int ri, int ci, int rf, int cf){
	if(cf == -1) cf=col-1;
	if(rf == -1) rf=row-1;
	Matrix _bm(rf-ri+1, cf-ci+1);
	for(int r=ri; r<=rf; r++)
		for(int c=ci; c<=cf; c++)
			_bm.set(r-ri, c-ci) = (*this)(r,c);
	return _bm;
}

void Matrix::setBlock(int ri, int ci, int rf, int cf, const Matrix& mblock){
	if(cf == -1) cf=col-1;
	if(rf == -1) rf=row-1;
	for(int r=ri; r<=rf; r++){
		for(int c=ci; c<=cf; c++){
			entry[r][c] = mblock.entry[r-ri][c-ci];
		}
	}
}

Matrix Matrix::getRow(int r){
	if(r == -1) r=row-1;
	return block(r,0,r,-1);
}
void Matrix::setRow(int r, const Matrix& mrow){
	if(r == -1) r=row-1;
	setBlock(r,0,r,-1,mrow);
}
void Matrix::setRow(int r, initializer_list< initializer_list<double> > iif){
	if(r == -1) r=row-1;
	Matrix miif = iif;
	setRow(r, miif);
}

// multiple rows
Matrix Matrix::getRow(int ri, int rf){
	if(rf == -1) rf=row-1;
	return block(ri,0,rf,-1);
}
void Matrix::setRow(int ri, int rf, const Matrix& mrow){
	if(rf == -1) rf=row-1;
	setBlock(ri,0,rf,-1,mrow);
}
void Matrix::setRow(int ri, int rf, initializer_list< initializer_list<double> > iif){
	if(rf == -1) rf=row-1;
	Matrix miif = iif;
	setRow(ri, rf, miif);
}

void Matrix::deleteRow(int ri, int rf){
	if(rf == -1) rf=col-1;
	Matrix storeTemp = (ri==0?
								(rf+1>col-1? 
											 Matrix(): getRow(rf+1,-1))
								:(rf+1>col-1?
											 getRow(0,ri-1):concat({
												{getRow(0,ri-1), getRow(rf+1,-1)}
											})
								)
						);

	*this = storeTemp;
}

void Matrix::deleteRow(int r){
	deleteRow(r,r);
}

// single coloumn
Matrix Matrix::getCol(int c){
	if(c == -1) c=col-1;
	return block(0,c,-1,c);
}
void Matrix::setCol(int c, const Matrix& mcol){
	if(c == -1) c=col-1;
	setBlock(0,c,-1,c,mcol);
}
void Matrix::setCol(int c, initializer_list < initializer_list<double> > iif){
	if(c == -1) c=col-1;
	Matrix miif = iif;
	setCol(c, miif);
}

// multiple coloumns
Matrix Matrix::getCol(int ci, int cf){
	if(cf == -1) cf=col-1;
	return block(0,ci,-1,cf);
}
void Matrix::setCol(int ci, int cf, const Matrix& mcol){
	if(cf == -1) cf=col-1;
	setBlock(0,ci,-1,cf,mcol);
}
void Matrix::setCol(int ci, int cf, initializer_list < initializer_list<double> > iif){
	if(cf == -1) cf=col-1;
	Matrix miif = iif;
	setCol(ci, cf, miif);
}

void Matrix::deleteCol(int ci, int cf){
	if(cf == -1) cf=col-1;
	Matrix storeTemp = (ci==0?
								(cf+1>col-1? 
											  Matrix(): getCol(cf+1,-1))
								:(cf+1>col-1?
											  getCol(0,ci-1):concat({
													{getCol(0,ci-1), getCol(cf+1,-1)}
											  })
								)
						);

	*this = storeTemp;
}

void Matrix::deleteCol(int c){
	deleteCol(c,c);
}

double Matrix::det(int method){
	if(row != col) return 0;
	
	double valdet(0.);
	switch(method){
		case 0:
			for(int c=0; c<col; c++){
				Matrix blkM;
				if(c == 0)
					blkM = this->block(1,1,row-1,col-1);
				else if(c == (col-1))
					blkM = this->block(1,0,row-1,col-2);
				else
					blkM = concat({{this->block(1,0,row-1,c-1), this->block(1,c+1,row-1,col-1)}});

				if(blkM.row == 1)
					valdet += pow(-1., (0+1)+(c+1)) * (*this)(0,c) * blkM(0,0);
				else
					valdet += pow(-1., (0+1)+(c+1)) * (*this)(0,c) * blkM.det();
			}
			break;
		case 1:
			std::unordered_map<std::string, Matrix> luofthis = LUdecomp();
			valdet = luofthis["U"].getDiag().mul_entry() * pow(-1., luofthis["row_change"](0));
			break;
	}

	return valdet;
}

std::unordered_map<std::string, Matrix> Matrix::LUdecomp(){
	Matrix PmatRight = Matrix::Identity(col);

	Matrix Umat = *this;

	if(!Umat.isZero()){
		for(int c=0; c<col; c++){
			Matrix PmatRightTemp = Matrix::Identity(col);
			if(Umat.getCol(c).isZero()){
				for(int cn=c+1; cn<col; cn++){
					if(!Umat.getCol(cn).isZero()){
						Matrix col_c = PmatRightTemp.getCol(c);
						Matrix col_cn = PmatRightTemp.getCol(cn);
						PmatRightTemp.setCol(c, col_cn);
						PmatRightTemp.setCol(cn, col_c);
						// PmatRightTemp.setCol(c, Matrix::Identity(col).getCol(cn));
						// PmatRightTemp.setCol(cn, Matrix::Identity(col).getCol(c));

						Umat = Umat * PmatRightTemp;
						PmatRight = PmatRight * PmatRightTemp;
						break;
					}
				}
			}

			if(Umat.getCol(c).isZero()) break;
		}
	}

	Matrix Pmat = Matrix::Identity(row);
	Matrix Lmat = Matrix::Identity(row);
	Matrix Emat = Matrix::Identity(row);

	int row_change = 0;

	Matrix Etemp = Matrix::Identity(row);

	bool isupper = isUpper();
	if(!isupper){

		while(true){
			int firstcol = 0;
			for(int c=firstcol; c<col-1; c++){
				for(int r=c+1; r<row; r++){
					Matrix Ptemp = Matrix::Identity(row);

					if(fabs(Umat(c,c)) <= EPSILON){
						for(int vr=c; vr<row; vr++){
							if(fabs(Umat(vr, c)) > EPSILON){
								Matrix row_vr = Ptemp.getRow(vr);
								Matrix row_c = Ptemp.getRow(c);

								Ptemp.setRow(c, row_vr);
								Ptemp.setRow(vr, row_c);
								
								Umat = Ptemp * Umat;
								Emat = Ptemp * Emat;
								Lmat = Lmat * Ptemp;
								Pmat = Ptemp * Pmat;

								row_change++;
								break;
							}				
						}
					}

					if(fabs(Umat(c,c)) <= EPSILON) break;

					Etemp.set(r,c) = double(-Umat(r,c) / Umat(c,c));
					Emat = Etemp * Emat;

					Umat = Etemp * Umat;

					// cout << "-------------------" << std::make_pair(r+1,c+1) << "---------------" << endl;
					// Etemp.print("Etemp");
					// Umat.print("Umat");

					Etemp.set(r,c) = -Etemp(r,c);
					Lmat = Lmat * Etemp;

					Etemp.set(r,c) = double(0.);
				}
			}

			n_rank_of_this_matrix = Umat.getDiag().numOfNonZero();

			if((n_rank_of_this_matrix < n_diag())
			   		&& !Umat.getRow(n_rank_of_this_matrix, -1).isZero()){
				for(int r=n_rank_of_this_matrix; r<row; r++){
					if(!Umat.getRow(r).isZero()){
						if(r == n_rank_of_this_matrix)
							break;
						else{
							Matrix Ptemp = Matrix::Identity(row);

							Matrix row_prev_last = Ptemp.getRow(n_rank_of_this_matrix);
							Matrix row_r = Ptemp.getRow(r);

							Ptemp.setRow(r, row_prev_last);
							Ptemp.setRow(n_rank_of_this_matrix, row_r);
							
							Umat = Ptemp * Umat;
							Emat = Ptemp * Emat;
							Lmat = Lmat * Ptemp;
							Pmat = Ptemp * Pmat;

							row_change++;
							break;
						}
					}
				}

				for(int c=0; c<col; c++){
					if(fabs(Umat(n_rank_of_this_matrix,c)) > EPSILON){
						Matrix PmatRightTemp = Matrix::Identity(col);

						Matrix col_prev_last = PmatRightTemp.getCol(n_rank_of_this_matrix);
						Matrix col_c = PmatRightTemp.getCol(c);

						PmatRightTemp.setCol(c, col_prev_last);
						PmatRightTemp.setCol(n_rank_of_this_matrix, col_c);

						Umat = Umat * PmatRightTemp;
						PmatRight = PmatRight * PmatRightTemp;

						row_change++;
						break;
					}
				}
				firstcol = n_rank_of_this_matrix;
			}else{
				Lmat = Pmat * Lmat;
				Emat = Emat * Pmat.t();
				break;
			}
		}
	}else{
		n_rank_of_this_matrix = Umat.getDiag().numOfNonZero();
	}

	// cout << (Pmat.t() * Lmat * Umat * PmatRight.t()).isSame(*this) << endl;

	std::unordered_map<std::string, Matrix> P_tLEU = {
		{"P", Pmat},
		{"PR", PmatRight}, // to change row of variable x in Ax -> M({P_R}x) -> Mz wit z = {P_R}x
		{"L", Lmat}, 
		{"E", Emat},
		{"U", Umat},
		{"row_change", row_change*Matrix::Ones(1,1)}
	};
	
	return P_tLEU;
}

std::unordered_map<std::string, Matrix> Matrix::rref(){
	std::unordered_map<std::string, Matrix> res_lu = LUdecomp();
	Matrix rref_mat = res_lu["U"];
	Matrix EPmat = res_lu["E"] * res_lu["P"];
	Matrix PRmat = res_lu["PR"];

	Matrix EuTemp = Matrix::Identity(EPmat.row);

	if(!isDiagonal()){
		for(int dc=n_rank_of_this_matrix-1; dc>0; dc--){
			for(int dr=dc-1; dr>=0; dr--){
				EuTemp.set(dr,dc) = double(-rref_mat(dr,dc)/rref_mat(dc,dc));

				EPmat = EuTemp * EPmat;
				rref_mat = EuTemp * rref_mat;

				cout << "-------------------" << std::make_pair(dr+1,dc+1) << "---------------" << endl;
				EuTemp.print("EuTemp");
				rref_mat.print("rref_mat");

				EuTemp.set(dr,dc) = double(0.);
			}
		}
	}

	Matrix diag_rref_mat = rref_mat.getDiag();

	Matrix invD;

	if(n_rank_of_this_matrix != 0){
		int ndiag = n_diag();
		if(n_rank_of_this_matrix < ndiag){
			invD = createDiag(concat({
							{1./diag_rref_mat.getRow(0, n_rank_of_this_matrix-1)},
							{Matrix::Ones(ndiag - n_rank_of_this_matrix, 1)}
						}));
		}else{
			Matrix _diag_part = 1./diag_rref_mat;
			if(row > col){
				invD = createDiag(concat({
					{_diag_part},
					{Matrix::Ones(row - n_rank_of_this_matrix, 1)}
				}));
			}else{
				invD = createDiag(_diag_part);
			}
		}
		
		EPmat = invD * EPmat;
		rref_mat = invD * rref_mat;
	}

	std::unordered_map<std::string, Matrix> res_rref = {
		{"E", EPmat},
		{"PR", PRmat},
		{"rref", rref_mat}
	};
	return res_rref;
}

std::unordered_map<std::string, Matrix> Matrix::null_space(){
	Matrix null_mat;
	std::unordered_map<std::string, Matrix> EnR = {
		{"E", Matrix::Identity(row)}
	};
	std::unordered_map<std::string, Matrix> res_ns;

	if(isZero())
		null_mat = Matrix::Identity(col);
	else{
		EnR = rref();

		if(n_rank_of_this_matrix == col)
			null_mat = Matrix::Zeros(col, 1);
		else{
			int n_depend = col-n_rank_of_this_matrix;

			Matrix I = Matrix::Identity(n_depend);

			Matrix F = EnR["rref"].block(0,n_rank_of_this_matrix,n_rank_of_this_matrix-1,-1);

			null_mat = concat({
				{-F},
				{I}
			});
			
			null_mat = EnR["PR"].t() * null_mat;
		}
	}

	res_ns = {
		{"E", EnR["E"]},
		{"ns", null_mat}
	};
	
	return res_ns;
}

std::unordered_map<std::string, Matrix> Matrix::qr(bool simple){
	std::unordered_map<std::string, Matrix> res_qr;
	Matrix Qort;

	if(simple){
		Qort = concat({{getCol(0).unit()}});
		for(int c=1; c<col; c++){
			Matrix vc = getCol(c);
			Matrix qn = getCol(c);
			for(int cq=0; cq<Qort.col; cq++){
				Matrix vcq = Qort.getCol(cq);
				qn = qn - dot(vcq, vc) * vcq;
			}
			Qort = concat({
				{Qort, qn.unit()}
			});
		}

		Matrix Rmat = Qort.t() * (*this);

		res_qr = {
			{"Q", Qort},
			{"R", Rmat}
		};		
	}else{
		std::unordered_map<std::string, Matrix> P_tLUPr_t = LUdecomp();
		
		Matrix P_tL = P_tLUPr_t["P"].t() * P_tLUPr_t["L"];

		// P_tLUPr_t["L"].print("P_tLUPr_t[L]");
		// P_tL.print("P_tL");

		Qort = concat({{P_tL.getCol(0).unit()}});
		for(int c=1; c<P_tL.col; c++){
			Matrix vc = P_tL.getCol(c);
			Matrix qn = P_tL.getCol(c);

			// qn.print("qn");

			for(int cq=0; cq<Qort.col; cq++){
				Matrix vcq = Qort.getCol(cq);
				qn = qn - dot(vcq, vc) * vcq;
			}
			Qort = concat({
				{Qort, qn.unit()}
			});
		}

		Matrix Pr = P_tLUPr_t["PR"];
		Matrix Pr_t = Pr.t();
		Matrix UPr_t = P_tLUPr_t["U"] * Pr_t;

		// cout << "in qr: " << (P_tL * UPr_t).isSame(*this) << endl;
		// (P_tL * UPr_t).print("------");

		Matrix Rfact = Qort.t() * P_tL;
		Matrix RofAll = Rfact * UPr_t;

		res_qr = {
			{"Q", Qort},
			{"R", RofAll},
			{"Pr", Pr}
		};
	}
	return res_qr;
}

std::unordered_map<std::string, Matrix> Matrix::qr_iter(double tol, int nloop, bool isprint){
	Matrix Ak = *this;
	Matrix Qk, Rk;
	Matrix Qtot = Matrix::Identity(Ak.row);
	while(true){
		std::unordered_map<std::string, Matrix> qrm = Ak.qr(true);
		Qk = qrm["Q"];
		Rk = qrm["R"];
		Ak = Rk * Qk;
		Qtot = Qtot * Qk;

		if(isprint)
			Ak.print("Ak");

		if(nloop != -1){
			if(nloop < 1) break;
			nloop--;
		}else if(Ak.isUpper(tol)) break;
	}

	std::unordered_map<std::string, Matrix> res_qr_iter = {
		{"Qtot", Qtot},
		{"U", Ak}
	};
	return res_qr_iter;
}

std::unordered_map<std::string, Matrix> Matrix::eig(bool compute_eigvec){
	Matrix eigval = qr_iter()["U"].getDiag();
	std::unordered_map<std::string, Matrix> res_eig = {
		{"eigval", eigval}
	};
	if(compute_eigvec){
		Matrix eigvec;
		Matrix detail(1,eigval.row);
		for(int r=0; r<eigval.row; r++){
			detail.set(r) = eigvec.col;
			eigvec = concat({{eigvec, ((*this) - eigval(r)*Matrix::Identity(row)).null_space()["ns"]}});
			detail.set(r) = eigvec.col - detail(r);
		}
		res_eig.emplace("eigvec", eigvec);
		res_eig.emplace("detail", detail);
	}
	return res_eig;
}

// std::unordered_map<std::string, Matrix> Matrix::SchurDecomp(){

// }

//////////////////////// OPERATOR
const Matrix& Matrix::operator=(const Matrix& m){ // Recommended
	row = m.row;
	col = m.col;

	entry = new double*[row];
	for(int r=0; r<row; r++){
		entry[r] = new double[col];
		for(int c=0; c<col; c++)
			entry[r][c] = m.entry[r][c];
	}
	return *this;
}

// void Matrix::operator=(const Matrix& m){
// 	row = m.row;
// 	col = m.col;
	
// 	entry = new double*[row];
// 	for(int r=0; r<row; r++){
// 		entry[r] = new double[col];
// 		for(int c=0; c<col; c++)
// 			entry[r][c] = m.entry[r][c];
// 	}
// }


Matrix Matrix::operator*(double a){
	Matrix mval(row, col);
	for(int r=0; r<row; r++)
		for(int c=0; c<col; c++)
			mval.set(r,c) = a * (*this)(r,c);
	return mval;
}

Matrix Matrix::operator*(const Matrix& M){
	Matrix resmul(row, M.col);
	for(int rs=0; rs<row; rs++){
		for(int cs=0; cs<M.col; cs++){
			for(int c=0; c<col; c++){
				resmul.set(rs,cs) += entry[rs][c]*M.entry[c][cs];
			}
		}
	}
	return resmul;
}

Matrix Matrix::operator+(const Matrix& M){
	Matrix mval(row, col);
	for(int r=0; r<row; r++)
		for(int c=0; c<col; c++)
			mval.set(r,c) = (*this)(r,c) + M.entry[r][c];
	return mval;
}

Matrix Matrix::operator-(const Matrix& M){
	Matrix mval(row, col);
	for(int r=0; r<row; r++)
		for(int c=0; c<col; c++)
			mval.set(r,c) = (*this)(r,c) - M.entry[r][c];
	return mval;
}

/////////////////////  STATIC
Matrix Matrix::MatrixNaN(int row_, int col_){
	Matrix mn(row_, col_);
	for(int r=0; r<row_; r++)
		for(int c=0; c<col_; c++)
			mn.set(r,c) = std::numeric_limits<double>::quiet_NaN();
	return mn;
}

Matrix Matrix::MatrixINF(int row_, int col_){
	Matrix mif(row_, col_);
	for(int r=0; r<row_; r++)
		for(int c=0; c<col_; c++)
			mif.set(r,c) = std::numeric_limits<double>::infinity();
	return mif;
}

Matrix Matrix::RandomUniformReal(int row_, int col_, Matrix bound){
	Matrix resrand(row_, col_);

	std::random_device rd;
	std::mt19937 gen(rd());

	for(int r=0; r<row_; r++){
		std::uniform_real_distribution<> dist(bound(r,0), bound(r,1));
		for(int c=0; c<col_; c++){
			resrand.set(r,c) = dist(gen);
		}
	}
	return resrand;
}

Matrix Matrix::Zeros(int row, int col){
	Matrix mstatic(row, col);
	for(int r=0; r<row; r++){
		for(int c=0; c<col; c++){
			mstatic.set(r,c) = double(0.);
		}
	}
	return mstatic;
}

Matrix Matrix::Zeros(size_m sm){
	return Matrix::Zeros(sm.first, sm.second);
}

Matrix Matrix::Ones(int row, int col){
	Matrix mstatic(row, col);
	for(int r=0; r<row; r++){
		for(int c=0; c<col; c++){
			mstatic.set(r,c) = 1.;
		}
	}
	return mstatic;
}

Matrix Matrix::Ones(size_m sm){
	return Matrix::Ones(sm.first, sm.second);
}

Matrix Matrix::Identity(int ndiag){
	Matrix mstatic = Matrix::Zeros(ndiag, ndiag);
	for(int d=0; d<ndiag; d++)
		mstatic.set(d,d) = 1.;
	return mstatic;
}

////////////////////////////////////////////////////  FRIEND FUNCTION
Matrix operator*(double scalar, Matrix M){
	return M * scalar;
}

Matrix operator/(double scalar, Matrix M){
	Matrix opm(M.size());
	for(int r=0; r<M.row; r++){
		for(int c=0; c<M.col; c++){
			opm.set(r,c) = 1./M(r,c);
		}
	}
	return opm;
}

////////////////////////////////////////////////////  OSTREAM
std::ostream& operator<<(std::ostream& os, const Matrix& M){
	for(int r=0; r< M.row; r++){
		for(int c=0; c< M.col; c++){
			os << M.entry[r][c] << ", ";
		}
		os << "\n";
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, std::pair<int, int> sm){
	os << "(" << sm.first << "," << sm.second << ")";
	return os;
}

std::ostream& operator<<(std::ostream& os, std::unordered_map<std::string, Matrix> um){
	typename std::unordered_map<std::string, Matrix>::iterator iter;
	for(iter = um.begin(); iter != um.end(); iter++){
		os << (iter->first) << ":\n" << round((iter->second), int(1./EPSILON)) << endl;
	}
	return os;
}

void Matrix::print(std::string _name_mat, bool print_round, int mult_round){
	if(print_round)
		cout << _name_mat << ":" << endl << round(*this, mult_round) << endl;
	else
		cout << _name_mat << ":" << endl << *this << endl;
}

// OTHER FUNCTION

double mul_entry(initializer_list< initializer_list <double> > iif){
	Matrix M = iif;
	return M.mul_entry();
}

double dot(Matrix v1, Matrix v2){
	Matrix v1dotv2 = v1.t() * v2;
	return v1dotv2(0,0);
}

Matrix createDiag(initializer_list< initializer_list <double> > iif){
	int col_ = (int)iif.begin()->size();
	Matrix MI = Matrix::Identity(col_);
	MI.setDiag(iif);
	return MI;
}

Matrix createDiag(const Matrix& mdiag){
	Matrix MI = Matrix::Identity(mdiag.row);
	MI.setDiag(mdiag);
	return MI;
}

Matrix round(Matrix M, int mult){
	Matrix _rounded(M.size());
	for(int r=0; r<M.row; r++){
		for(int c=0; c<M.col; c++){
			if(std::isnan(M(r,c)))
				_rounded.set(r,c) = M(r,c);
			else
				_rounded.set(r,c) = int(M(r,c)*mult)/static_cast<double>(mult);
		}
	}
	return _rounded;
}

Matrix concat(const initializer_list< initializer_list <Matrix> >& iim){
	int nrowcon = (int)iim.size();
	int ncolcon = (int)iim.begin()->size();
	
	int rown(0), coln(0);
	for(int rc=0; rc<nrowcon; rc++){
		int rowpart;
		for(int cc=0; cc<ncolcon; cc++){
			rowpart = (iim.begin() + rc)->begin()[cc].row;
			if(rowpart != 0) break;
		}
		rown += rowpart;
	}

	for(int cc=0; cc<ncolcon; cc++){
		int colpart;
		for(int rr=0; rr<nrowcon; rr++){
			colpart = (iim.begin() + 0)->begin()[cc].col;
			if(colpart != 0) break;
		}
		coln += colpart;
	}

	Matrix rescon(rown, coln);

	for(int rc=0, cor=0; rc<nrowcon; rc++){
		int row_;
		for(int cc=0; cc<ncolcon; cc++){
			row_ = (int)(iim.begin() + rc)->begin()[cc].row;
			if(row_ != 0) break;
		}

		ncolcon = (int)(iim.begin() + rc)->size();
		for(int cc=0, coc=0; cc<ncolcon; cc++){
			int col_ = (int)(iim.begin() + rc)->begin()[cc].col;
			if(col_ == 0) continue;
			rescon.setBlock(cor, coc, cor+row_-1, coc+col_-1, (iim.begin() + rc)->begin()[cc]);
			coc += col_;
		}

		cor += row_;
	}
	return rescon;
}

std::unordered_map<std::string, Matrix > solve(Matrix A, Matrix b){
	std::unordered_map<std::string, Matrix> res_ns_A = A.null_space();
	Matrix Erref = res_ns_A["E"];
	Matrix rref_b = Erref * b;

	Matrix ns_A = res_ns_A["ns"];

	Matrix ps_A = rref_b;
	std::unordered_map<std::string, Matrix > res_solve;

	for(int c=0; c<b.col; c++){
		if(rref_b.getCol(c).numOfNonZero() > A.n_rank_of_this_matrix){
			ps_A.setCol(c, Matrix::MatrixNaN(b.row,1));
			cout << "\nWARNING: coloumn-" << c+1 << " of matrix b, is not inside image(A), thus no solution" << endl;
		}
	}

	res_solve = {
		{"ps", ps_A},
		{"ns", ns_A}
	};
	return res_solve;
}