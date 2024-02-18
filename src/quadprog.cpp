#include "mex.hpp"
#include "mexAdapter.hpp"

#include "utilities.hpp"
#include "sparse.hpp"
#include "Highs.h"
#include "limits"

class MexFunction
    : public matlab::mex::Function
{
public:
    MexFunction() { matlabPtr = getEngine(); }
    ~MexFunction() {}
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        
        if (inputs.size() < 4)
            utilities::error("At least 4 inputs are required: Q, f, A, b.");

        matlab::data::Array Qmat = std::move(inputs[0]);
        matlab::data::TypedArray<double> f = std::move(inputs[1]);
        matlab::data::Array Amat = std::move(inputs[2]);
        matlab::data::TypedArray<double> b = std::move(inputs[3]);
        
        if (Amat.getDimensions().size() != 2)
            utilities::error("A must be a Matrix.");

        std::size_t nx = f.getNumberOfElements();
        std::size_t nc = Amat.getDimensions()[0];
        
        if (b.getDimensions().size() != 2 || b.getDimensions()[0] != nc || b.getDimensions()[1] != 2)
            utilities::error("b must be a {} x 2 array with [lb,ub] for each constraint.",nc);

        if (Amat.getDimensions()[1] != nx)
            utilities::error("The number of columns in A must be equal to the number of elements in f.");
        
        utilities::Sparse<double> A(nc, nx);
        utilities::Sparse<double> Q(nx, nx);
        A.set(Amat);
        Q.set(Qmat);
        
        HighsModel model;
        
        model.hessian_.dim_ = nx;
        model.hessian_.format_ = HessianFormat::kTriangular;
        model.hessian_.start_.resize(nx + 1);
        model.hessian_.index_.resize(Q.getNumberOfNonZeroElements());
        model.hessian_.value_.resize(Q.getNumberOfNonZeroElements());
        Q.getCsc<HighsInt>(model.hessian_.start_, model.hessian_.index_, model.hessian_.value_);

        model.lp_.num_col_ = nx;
        model.lp_.num_row_ = nc;
        model.lp_.sense_ = ObjSense::kMinimize;
        model.lp_.offset_ = 0.;

        model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
        model.lp_.a_matrix_.start_.resize(nx + 1);
        model.lp_.a_matrix_.index_.resize(A.getNumberOfNonZeroElements());
        model.lp_.a_matrix_.value_.resize(A.getNumberOfNonZeroElements());

        A.getCsc<HighsInt>(model.lp_.a_matrix_.start_, model.lp_.a_matrix_.index_, model.lp_.a_matrix_.value_);

        model.lp_.col_cost_.resize(nx);
        std::copy(f.cbegin(), f.cend(), model.lp_.col_cost_.begin());

        model.lp_.row_lower_.resize(nc);
        model.lp_.row_upper_.resize(nc);
        for (std::size_t iRow = 0; iRow < nc; iRow++)
        {
            model.lp_.row_lower_[iRow] = b[iRow][0];
            model.lp_.row_upper_[iRow] = b[iRow][1];
        }

        model.lp_.col_lower_.resize(nx);
        model.lp_.col_upper_.resize(nx);
        if (inputs.size() > 4)
        {
            matlab::data::TypedArray<double> bnds = std::move(inputs[4]);
            if (bnds.getDimensions().size() != 2 || bnds.getDimensions()[0] != nx)
                utilities::error("bnds must be a {} x 2 array with [lb,ub] for each variable.",nx);
            for (std::size_t iVar = 0; iVar < nx; iVar++) {
                model.lp_.col_lower_[iVar] = bnds[iVar][0];
                model.lp_.col_upper_[iVar] = bnds[iVar][1];
            }
        }
        else {
            std::fill(model.lp_.col_lower_.begin(), model.lp_.col_lower_.end(), -kHighsInf);
            std::fill(model.lp_.col_upper_.begin(), model.lp_.col_upper_.end(), kHighsInf);
        }

        Highs highs;
        HighsStatus status = highs.passModel(model);
        if (status != HighsStatus::kOk)
            utilities::error("Failed to pass the model to Highs.");
        
        status = highs.run();
        if (status != HighsStatus::kOk)
            utilities::error("Failed to solve the problem.");
        
        const HighsModelStatus &model_status = highs.getModelStatus();
        const HighsInfo &info = highs.getInfo();
        const HighsSolution &solution = highs.getSolution();
        const HighsBasis& basis = highs.getBasis();

        matlab::data::ArrayFactory factory;
        matlab::data::TypedArray<double> x = factory.createArray<double>({nx, 1});
        std::copy(solution.col_value.begin(), solution.col_value.end(), x.begin());

        outputs[0] = std::move(x);
        if (outputs.size() > 1) {
            matlab::data::TypedArray<double> fVal = factory.createScalar(info.objective_function_value);
            outputs[1] = std::move(fVal);
        }

        if (outputs.size() > 2) {
            matlab::data::StructArray infoStruct = factory.createStructArray({1, 1}, {"column_dual", "column_status", "row_dual", "row_status", "row_value", "status"});
            matlab::data::TypedArray<double> column_dual = factory.createArray<double>({nx, 1});
            matlab::data::TypedArray<matlab::data::MATLABString> column_status = factory.createArray<matlab::data::MATLABString>({nx, 1});
            matlab::data::TypedArray<double> row_dual = factory.createArray<double>({nc, 1});
            matlab::data::TypedArray<double> row_value = factory.createArray<double>({nc, 1});
            matlab::data::TypedArray<matlab::data::MATLABString> row_status = factory.createArray<matlab::data::MATLABString>({nc, 1});
            matlab::data::TypedArray<matlab::data::MATLABString> status = factory.createScalar(highs.modelStatusToString(model_status));
            std::copy(solution.col_dual.begin(), solution.col_dual.end(), column_dual.begin());
            std::copy(solution.row_value.cbegin(), solution.row_value.cend(), row_value.begin());
            std::copy(solution.row_dual.cbegin(), solution.row_dual.cend(), row_dual.begin());
            std::transform(basis.col_status.cbegin(), basis.col_status.cend(), column_status.begin(), [&highs](const HighsBasisStatus &status) { return highs.basisStatusToString(status); });
            std::transform(basis.row_status.cbegin(), basis.row_status.cend(), row_status.begin(), [&highs](const HighsBasisStatus &status) { return highs.basisStatusToString(status); });

            infoStruct[0]["column_dual"] = std::move(column_dual);
            infoStruct[0]["column_status"] = std::move(column_status);
            infoStruct[0]["row_dual"] = std::move(row_dual);
            infoStruct[0]["row_status"] = std::move(row_status);
            infoStruct[0]["row_value"] = std::move(row_value);
            infoStruct[0]["status"] = std::move(status);
            outputs[2] = std::move(infoStruct);
        }
    }
};