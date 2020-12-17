# pragma once

namespace scd {
    namespace test {
// Code based on
// https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix

        template<typename M>
        M ReadCsvMatrix(const std::string &path, uint line_start) {
            std::ifstream indata;
            indata.open(path);
            if (indata.fail())
                std::cerr << "Failed to open file at " << path << std::endl;
            std::string line;
            std::vector<double> values;
            uint line_n = 0;
            uint line_end = line_start + M::RowsAtCompileTime - 1;
            while ((std::getline(indata, line)) && (line_n <= line_end)) {
                if (line_n >= line_start) {
                    std::stringstream lineStream(line);
                    std::string cell;
                    while (std::getline(lineStream, cell, ',')) {
                        values.push_back(std::stod(cell));
                    }
                }
                ++line_n;
            }

            return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime,
                    Eigen::RowMajor>>(values.data(), M::RowsAtCompileTime, M::ColsAtCompileTime);
        }

// TODO (jay) Combine with read_csv_matrix(). Currently separate because RowMajor not allowed for column vectors.
        template<typename M>
        M ReadCsvVector(const std::string &path, uint line_start) {
            std::ifstream indata;
            indata.open(path);
            if (indata.fail())
                std::cerr << "Failed to open file at " << path << std::endl;
            std::string line;
            std::vector<double> values;
            uint line_n = 0;
            uint line_end = line_start + M::RowsAtCompileTime - 1;
            while ((std::getline(indata, line)) && (line_n <= line_end)) {
                if (line_n >= line_start) {
                    std::stringstream lineStream(line);
                    std::string cell;
                    while (std::getline(lineStream, cell, ',')) {
                        values.push_back(std::stod(cell));
                    }
                }
                ++line_n;
            }

            return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime,
                    Eigen::ColMajor>>(values.data(), M::RowsAtCompileTime, M::ColsAtCompileTime);
        }

// TODO (jay) Refactor this using templates instead of overloading
        void ReadShapeFile(std::string file_path, uint n_line, Superquadric<3>& SQ){
            const uint n = 3;
            std::ifstream indata;
            indata.open(file_path);

            Eigen::Matrix<double,1,n> a = scd::test::ReadCsvMatrix<Eigen::Matrix<double,1,n>>(file_path, n_line);
            Eigen::Matrix<double,1,n-1> e = scd::test::ReadCsvMatrix<Eigen::Matrix<double,1,n-1>>(file_path, n_line+1);
            Eigen::Matrix<double,n,n> R = scd::test::ReadCsvMatrix<Eigen::Matrix<double,n,n>>(file_path, n_line+2);
            Eigen::Matrix<double,n,1> t = scd::test::ReadCsvVector<Eigen::Matrix<double,n,1>>(file_path, n_line+5);

            SQ.a = {a[0], a[1], a[2]};
            SQ.e = {e[0], e[1]};

            Eigen::Matrix<double,4,4> X;
            X.setIdentity();
            X.block<3,3>(0,0) = R;
            X.block<3,1>(0,3) = t;
            SQ.X.matrix() =X;
            return;
        }

        void ReadTestFile(std::string file_path, uint n_line, Superquadric<3>& SQ1, Superquadric<3>& E2, Superquadric<3>& E2c, CollideResult<3>& result){
            std::ifstream indata;
            indata.open(file_path);
            if (indata.fail())
                std::cerr << "Failed to open file " << file_path << std::endl;

            const uint n = 3;
            ReadShapeFile(file_path, n_line, SQ1);
            ReadShapeFile(file_path, n_line+8, E2);
            ReadShapeFile(file_path, n_line+16, E2c);
            Eigen::Matrix<double,1,3> d = scd::test::ReadCsvMatrix<Eigen::Matrix<double,1,3>>(file_path, n_line+24);
            result.angles[0] = d(0);
            result.angles[1] = d(1);
            result.collide = static_cast<bool>(d(2));
            return;
        }

    } // namespace test
} // namespace scd