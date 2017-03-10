#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

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
	virtual TransferMatrix getTransferMatrixAtFrequency(double f) = 0;
};

class Resistor : public TwoPortElement {
public:
	Resistor(double r) : resistance(r) {}

	virtual TransferMatrix getTransferMatrixAtFrequency(double f) = 0;

	double resistance;
};

class SeriesResistor : public Resistor {
public:
	SeriesResistor(double r) : Resistor(r) {}

	TransferMatrix getTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, -resistance, 0, 1);
	}
};

class ShuntResistor : public Resistor {
public:
	ShuntResistor(double r) : Resistor(r) {}

	TransferMatrix getTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, 0, -1.0/resistance, 1);
	}
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

	virtual TransferMatrix getTransferMatrixAtFrequency(double f) = 0;

	double capacitance; // in pF
};

class SeriesCapacitor : public Capacitor {
public:
	SeriesCapacitor(double c) : Capacitor(c) {}

	TransferMatrix getTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, ComplexNumber(0, -xc(capacitance, f)), 0, 1);
	}
};

class ShuntCapacitor : public Capacitor {
public:
	ShuntCapacitor(double c) : Capacitor(c) {}

	TransferMatrix getTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, 0, ComplexNumber(0, -1.0/xc(capacitance, f)), 1);
	}
};



class Inductor : public TwoPortElement {
public:
	Inductor(double l) : inductance(l) {}

	virtual TransferMatrix getTransferMatrixAtFrequency(double f) = 0;

	double inductance; // in nanohenries
};

class SeriesInductor : public Inductor {
public:
	SeriesInductor(double l) : Inductor(l) {}

	TransferMatrix getTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, ComplexNumber(0, -xl(inductance/1000.0, f)), 0, 1);
	}
};

class ShuntInductor : public Inductor {
public:
	ShuntInductor(double l) : Inductor(l) {}

	TransferMatrix getTransferMatrixAtFrequency(double f) override {
		return TransferMatrix(1, 0, ComplexNumber(0, -1.0/xl(inductance/1000.0, f)), 1);
	}
};

class TwoPortNetwork : public TwoPortElement {
public:
	void add(TwoPortElement* element) {
		elements.push_back(element);
	}

	TransferMatrix getTransferMatrixAtFrequency(double f) override {
		TransferMatrix result(1, 0, 0, 1); // identity complex matrix

		if (!elements.empty()) {
			for (auto& e : elements) {
				result = result * e->getTransferMatrixAtFrequency(f);
			}
		}

		return result;
	}

	std::vector<TwoPortElement*> elements;
};


int main() {

	TwoPortNetwork network;
	network.add(new SeriesResistor(50)); // source R
	network.add(new ShuntCapacitor(80.35));
	network.add(new SeriesInductor(44.11));
	network.add(new ShuntCapacitor(106.16));
	network.add(new SeriesInductor(46.15));
	network.add(new ShuntCapacitor(107.48));
	network.add(new SeriesInductor(46.15));
	network.add(new ShuntCapacitor(106.16));
	network.add(new SeriesInductor(44.11));
	network.add(new ShuntCapacitor(80.35));
	network.add(new ShuntResistor(50)); // load R

	for (double frequency = 0; frequency <= 220000000; frequency += 200000) {
		TransferMatrix filter = network.getTransferMatrixAtFrequency(frequency);

		double outVoltage = 1.0/magnitude(filter.A);
		double outGain = 20.0 * log10(outVoltage); // out/in
		double outPhase = phase(filter.A);

		std::cout << (frequency/1000000.0) << "\t" << outGain << "\t" << outPhase << std::endl;
	}
}