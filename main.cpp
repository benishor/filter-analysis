#include <iostream>
#include <iomanip>
#include <cmath>

class ComplexNumber {
public:
	ComplexNumber(double r)
		: real(r), imag(0.0) {}

	ComplexNumber(double r, double i)
		: real(r), imag(i) {}


	double real, imag;
};

std::ostream& operator << (std::ostream& out, const ComplexNumber& c) {
	out << std::setprecision(6) << "(" << c.real << (c.imag < 0 ? "" : "+") << c.imag << "i" << ")";
	return out;
}

ComplexNumber operator * (ComplexNumber a, ComplexNumber b) {
	return ComplexNumber(
		a.real*b.real - a.imag*b.imag,
		a.real*b.imag + a.imag*b.real
	);
}

ComplexNumber operator + (ComplexNumber a, ComplexNumber b) {
	return ComplexNumber(a.real+b.real, a.imag+b.imag);
}

double magnitude(ComplexNumber a) {
	return sqrt(a.real*a.real + a.imag*a.imag);
}

double phase(ComplexNumber a) {
	double angleInRadians = atan2(a.imag, a.real);
	return (angleInRadians > 0 ? angleInRadians : (2.0*M_PI + angleInRadians)) * 360.0 / (2.0 * M_PI);
}

class VoltageAndCurrent {
public:
	VoltageAndCurrent(ComplexNumber voltage, ComplexNumber current)
		: v(voltage), i(current) {}

	ComplexNumber v;
	ComplexNumber i;
};


class TransferMatrix {
public:
	TransferMatrix() 
		: A(1), B(0), C(0), D(1) {}
	TransferMatrix(ComplexNumber a, ComplexNumber b, ComplexNumber c, ComplexNumber d)
		: A(a), B(b), C(c), D(d) {}

	ComplexNumber A, B, C, D;
};

std::ostream& operator << (std::ostream& out, const TransferMatrix& m) {
	out << "[" << m.A << ", " << m.B << "]" << std::endl;
	out << "[" << m.C << ", " << m.D << "]" << std::endl;
	return out;
}

TransferMatrix operator * (TransferMatrix a, TransferMatrix b) {
	return TransferMatrix(
		a.A*b.A+a.B*b.C, a.A*b.B+a.B*b.D,
		a.C*b.A+a.D*b.C, a.C*b.B+a.D*b.D
	);
}

VoltageAndCurrent operator * (TransferMatrix matrix, VoltageAndCurrent input) {
	return VoltageAndCurrent(
		matrix.A*input.v + matrix.B*input.i,
		matrix.C*input.v + matrix.D*input.i
	);
}

std::ostream& operator << (std::ostream& out, const VoltageAndCurrent& vi) {
	out << "[" << vi.v << "]" << std::endl;
	out << "[" << vi.i << "]" << std::endl;
	return out;
}


class TwoPortElement {
public:
	virtual TransferMatrix getSeriesTransferMatrixAtFrequency(double f) = 0;
	virtual TransferMatrix getShuntTransferMatrixAtFrequency(double f) = 0;
};

class Resistor : public TwoPortElement {
public:

	Resistor(double r) : resistance(r) {}

	TransferMatrix getSeriesTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, -resistance, 0, 1);
	}

	TransferMatrix getShuntTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, 0, -1.0/resistance, 1);
	}

	double resistance;
};

double xc(double picofarads, double hz) {
	return 1.0/(2.0 * M_PI * hz * picofarads/1000000000000.0);
}

double xl(double microhenries, double hz) {
	return 2.0 * M_PI * hz * microhenries/1000000.0;
}

class Capacitor : public TwoPortElement {
public:

	Capacitor(double c) : capacitance(c) {}

	TransferMatrix getSeriesTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, ComplexNumber(0, -xc(capacitance, f)), 0, 1);
	}

	TransferMatrix getShuntTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, 0, ComplexNumber(0, -1.0/xc(capacitance, f)), 1);
	}

	double capacitance; // in pF
};

class Inductor : public TwoPortElement {
public:

	Inductor(double l) : inductance(l) {}

	TransferMatrix getSeriesTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, ComplexNumber(0, -xl(inductance/1000.0, f)), 0, 1);
	}

	TransferMatrix getShuntTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, 0, ComplexNumber(0, -1.0/xl(inductance/1000.0, f)), 1);
	}

	double inductance; // in nanohenries
};


int main() {
	// ComplexNumber c(1, 2.2);
	// std::cout << "Hello, world! Here's a complex number: " << c << std::endl;

	// Resistor r(100);
	// TransferMatrix m1 = r.getSeriesTransferMatrixAtFrequency(100);
	// TransferMatrix m2 = r.getShuntTransferMatrixAtFrequency(100);
	// std::cout << "Here's a series transfer matrix: " << std::endl << m1 << std::endl;
	// std::cout << "Here's a shunt transfer matrix: " << std::endl << m2 << std::endl;

	// Capacitor cap(100);
	// TransferMatrix m3 = cap.getSeriesTransferMatrixAtFrequency(1000000);
	// TransferMatrix m4 = cap.getShuntTransferMatrixAtFrequency(1000000);

	// std::cout << "Here's a series capacitor transfer matrix: " << std::endl << m3 << std::endl;
	// std::cout << "Here's a shunt capacitor transfer matrix: " << std::endl << m4 << std::endl;

	// double capacitance = 100.0;
	// std::cout << "Xc: " << xc(frequency, capacitance) << std::endl;

	// TransferMatrix a(1, 2, 3, 4), b(5, 6, 7, 8);
	// TransferMatrix cc = a*b;
	// std::cout << "Multiplication result: " << std::endl << cc << std::endl;

	// ----

	Capacitor c1(80.35);
	Inductor l1(44.11);
	Capacitor c2(106.16);
	Inductor l2(46.15);
	Capacitor c3(107.48);
	Inductor l3(46.15);
	Capacitor c4(106.16);
	Inductor l4(44.11);
	Capacitor c5(80.35);
	Resistor r1(50);

	for (double frequency = 0; frequency <= 220000000; frequency += 200000) {
		TransferMatrix filter = r1.getSeriesTransferMatrixAtFrequency(frequency) *
								c1.getShuntTransferMatrixAtFrequency(frequency) * 
								l1.getSeriesTransferMatrixAtFrequency(frequency) *
								c2.getShuntTransferMatrixAtFrequency(frequency) * 
								l2.getSeriesTransferMatrixAtFrequency(frequency) *
								c3.getShuntTransferMatrixAtFrequency(frequency) * 
								l3.getSeriesTransferMatrixAtFrequency(frequency) *
								c4.getShuntTransferMatrixAtFrequency(frequency) * 
								l4.getSeriesTransferMatrixAtFrequency(frequency) *
								c5.getShuntTransferMatrixAtFrequency(frequency) * 
								r1.getShuntTransferMatrixAtFrequency(frequency);

		// std::cout << "A is: " << std::endl << tm1 << std::endl;
		// std::cout << "B is: " << std::endl << tm2 << std::endl;

		// std::cout << "Filter matrix is: " << std::endl << filter << std::endl;
		// double output = 1.0 / magnitude(filter.A);
		// std::cout << frequency << "\t"  << std::fixed << std::setprecision(6) <<  output << std::setprecision(0) << std::endl;

		// VoltageAndCurrent input(1, 0);
		// VoltageAndCurrent output = filter * input;


		double outVoltage = 1.0/magnitude(filter.A);
		double outGain = 20.0 * log10(outVoltage); // out/in
		double outPhase = phase(filter.A);
		// std::cout << "Filter output at " << frequency << " Hz is: " << std::fixed << std::setprecision(6) << outVoltage << " (" << outGain << " dB)" << std::endl;
		// std::cout << "Filter output phase at " << frequency << " Hz is: " << std::fixed << std::setprecision(6) << phase(filter.A) << std::endl;
		std::cout << (frequency/1000000.0) << "\t" << outGain << "\t" << outPhase << std::endl;
	}

	// std::cout << "Xc= " << xc(47000, 10000) << std::endl;

	// Resistor r1(100), r2(200);
	// for (frequency = 100; frequency < 1000; frequency += 100) {
	// 	TransferMatrix filter = r1.getSeriesTransferMatrixAtFrequency(frequency) * 
	// 							r2.getShuntTransferMatrixAtFrequency(frequency);
	// 	// std::cout << "Filter matrix is: " << std::endl << filter << std::endl;
	// 	double output = 1.0 / magnitude(filter.A);
	// 	std::cout << frequency << "\t"  << std::fixed << std::setprecision(6) <<  output << std::setprecision(0) << std::endl;
	// }

	// std::cout << "Series transfer matrix is: " << std::endl << rc1.getSeriesTransferMatrixAtFrequency(frequency) << std::endl;
	// std::cout << "Shunt transfer matrix is: " << std::endl << rc1.getShuntTransferMatrixAtFrequency(frequency) << std::endl;
	// std::cout << "Filter matrix is: " << std::endl << filter << std::endl;

	// VoltageAndCurrent input(1, 0);
	// VoltageAndCurrent output = filter * input;

	// std::cout << "Filter output at " << frequency << " Hz is: " << std::endl << output << std::endl;


	// TransferMatrix n12(
	// 	ComplexNumber(-7,0), ComplexNumber(0, -4), 
	// 	ComplexNumber(0, 2), ComplexNumber(1, 0)
	// );
	// TransferMatrix n3(
	// 	ComplexNumber(1,0), ComplexNumber(0, 3), 
	// 	ComplexNumber(0, 0), ComplexNumber(1, 0)
	// );

	// TransferMatrix n4(
	// 	ComplexNumber(1,0), ComplexNumber(0, 0), 
	// 	ComplexNumber(0, -2), ComplexNumber(1, 0)
	// );

	// TransferMatrix n123 = n12*n3;
	// std::cout << "n123 matrix is: " << std::endl << n123 << std::endl;

	// TransferMatrix n1234 = n123*n4;
	// std::cout << "n1234 matrix is: " << std::endl << n1234 << std::endl;

	// ComplexNumber cc1(0, 2), cc2(0, 3);
	// std::cout << "Multiplication result is : " << std::endl << cc1*cc2 << std::endl;

	// ComplexNumber cc3(1, 0), cc4(1, 0);
	// std::cout << "Multiplication result is : " << std::endl << cc3*cc4 << std::endl;

	// ComplexNumber sum = cc1*cc2 + cc3*cc4;
	// std::cout << "Result is : " << std::endl << sum << std::endl;

}