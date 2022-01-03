#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

//����� ��� ��������� ��������� �����
class RandomGenerator
{
public:
	static std::mt19937 & getMt19937();

private:
	RandomGenerator();
	~RandomGenerator() {}
	static RandomGenerator& instance();

	RandomGenerator(RandomGenerator const&) = delete;
	RandomGenerator& operator= (RandomGenerator const&) = delete;

	std::mt19937 mMt;
};

RandomGenerator::RandomGenerator() {
	std::random_device rd;

	if (rd.entropy() != 0) {
		std::seed_seq seed{ rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
		mMt.seed(seed);
	}
	else {
		auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		mMt.seed(seed);
	}
}

RandomGenerator& RandomGenerator::instance() {
	static RandomGenerator s;
	return s;
}

std::mt19937 & RandomGenerator::getMt19937() {
	return RandomGenerator::instance().mMt;
}

//��� �������

//��������� �������
vector <double> sampleGeneration(int n)
{
    std::mt19937 &mt = RandomGenerator::getMt19937();
	std::uniform_real_distribution<double> dist(0.0, 1.0);
    vector <double> sample;
    sample.resize(n);
    int i = 0;
    while (i < n)
    {
        sample[i] = dist(mt);
        i++;
    }

    return sample;
}

//������������� ������� � ������� ������������� �����
vector <double> modelRayleighRaspr(int n)
{
    vector <double> u1 = sampleGeneration(n);
    vector <double> u2 = sampleGeneration(n);
    vector <double> u3 = sampleGeneration(n);
    vector <double> u4 = sampleGeneration(n);
    vector <double> x;
    x.resize(u1.size());
    vector <double> y;
    y.resize(u3.size());

    for (int i = 0; i < x.size(); i++)
    {
        x[i] = pow(-2.0*log(u1[i]), 0.5) * cos(2*3.142*u2[i]);
        y[i] = pow(-2.0*log(u3[i]), 0.5) * cos(2*3.142*u4[i]);
    }

    vector <double> z;
    z.resize(x.size());

    for (int i = 0; i < z.size(); i++)
    {
        z[i] = pow(pow(x[i], 2.0) + pow(y[i], 2.0), 0.5);
    }

    return z;
}

//������������� ����������������� �������������
vector <double> modelExpRaspr(int n, double a = 1.0)
{
    vector <double> u = sampleGeneration(n);
    vector <double> x;
    x.resize(u.size());

    for (int i = 0; i < x.size(); i++)
    {
        x[i] = -(1/a) * log(1 - u[i]);
    }

    return x;
}

//������� ��������� ����������������� �������������
double expRaspr(double x, double m = 1.0, double s = 0)
{
    if(x >= 0)
        return (1 - exp(-m*x));
    else return 0;
}

//������� ��������� ������������� �����
double relRaspr(double x, double m = 1.0, double s = 0)
{
    if(x >= 0)
        return (1 - exp(-(pow(x, 2) / (2 * pow(m, 2)))));
    else return 0;
}

//�������� ��-������� �������
double pirsonSquare(vector <double> &sample, int k, char r)
{
    //���������� ������� ��� �������������
    sort(sample.begin(), sample.begin()+sample.size());

    //��� 1: ������������� ������
    //��� 1.1: ������ �������
    double max = *max_element(sample.begin(), sample.end());
    double min = *min_element(sample.begin(), sample.end());
    double d = max - min;

    double h = d / k;


    //��� 1.2: ��������� �� ���������
    //������ ��� �������� ��������� � ������ ��������
    vector <double> kLength;
    kLength.resize(k);

    for(int i = 0, j = 0; i < sample.size()-1; i++)
    {
        if(sample[i] <= (sample[0] + h))
        {
            //���� ��������� �������� �������� � ��������, �� ������� ��������� � ������ �������� ������������� �� 1
            kLength[j]++;

            //���������� ���� �� ���������� ��������
            h = d / k;
            //����� j �� ������� ���������
            j = 0;
        }

        else
        {
            //���� ��������� �������� �� �������� � ��������, �� �������� �����������
            i--;
            //���������� j �� ���������� ���������
            j++;
            //���������� h �� ���������� ���������
            h += d / k;
        }
    }

    //����������� ���� � ��������� ��������� � ����� ����
    kLength[k-1]++;

    //���������� �������� ����
    h = d / k;

    //��� 1.3: ����������� ��������� � ��������
    //������ ��� �������� �����������
    vector <double> pInterval;
    pInterval.resize(k);

    for(int i = 0; i < k; i++)
    {
        pInterval[i] = kLength[i] / sample.size();
    }

    //������ ��� �������� �������� ��������� �����
    vector <double> xInterval;
    xInterval.resize(k);

    //cout << "Boundary points: " << endl;
    xInterval[0] = sample[0] + h;
    for(int i = 1; i < k-1; i++)
    {
        xInterval[i] = xInterval[i-1] + h;
        //if(i == 1)
            //cout << xInterval[i-1] << endl;
        //cout << xInterval[i] << endl;
    }

    //cout << endl;

    //��� 2: ��-������� �������
    //��� 2.1: ���������� ��������� ������� ��������� � ��������
    //������ ��� �������� ����������� ��������� � ��������
    vector <double> pT;
    pT.resize(k);

    double sum = 0;
    //cout << "Expected frequency of hitting the interval: " << endl;
    switch(r)
    {
    case 'e':
        for(int i = 0; i < k; i++)
        {
            if(i == 0)
                pT[i] = expRaspr(xInterval[i]);
            else if(i > 0 && i < k - 1)
                pT[i] = expRaspr(xInterval[i]) - sum;
            else pT[i] = 1 - sum;
            //cout << pT[i] << endl;
            sum += pT[i];
        }
        break;
    case 'r':
        for(int i = 0; i < k; i++)
        {
            if(i == 0)
                pT[i] = relRaspr(xInterval[i]);
            else if(i > 0 && i < k - 1)
                pT[i] = relRaspr(xInterval[i]) - sum;
            else pT[i] = 1 - sum;
            //cout << pT[i] << endl;
            sum += pT[i];
        }
        break;
    }

    //cout << endl;

    //��� 2.2: ���������� ����������
    double pSquare = 0;
    for(int i = 0; i < k; i++)
    {
        pSquare += pow((pInterval[i] - pT[i]), 2) / pT[i];
    }

    return sample.size() * pSquare;
}

double monteKarlo(double xPSquare, int n, int k, char r, int N = 16600)
{
    double m = 0;
    double yPSquare;
    cout << "Waiting...\n\n";
    switch(r)
    {
        case 'e':
            for(int i = 0; i < N; i++)
            {
                vector <double> y = modelExpRaspr(n);
                yPSquare = pirsonSquare(y, k, r);
                if(yPSquare > xPSquare)
                m = m + 1;
            }
            break;

        case 'r':
            for(int i = 0; i < N; i++)
            {
                vector <double> y = modelRayleighRaspr(n);
                yPSquare = pirsonSquare(y, k, r);
                if(yPSquare > xPSquare)
                m = m + 1;
            }

    }

    return m/N;
}

void pValueTab(double pSquare, int k, double aValue = 0.05)
{
    vector <double> pValue {0.01, 0.025, 0.05, 0.95, 0.975, 0.99};
    vector <double> pStatic7 {18.5, 16.0, 14.1, 2.17, 1.69, 1.24};
    vector <double> pStatic9 {21.7, 19.0, 16.9, 3.33, 2.7, 2.09};
    vector <double> pStatic13 {27.7, 24.7, 22.4, 5.89, 5.01, 4.11};

    double pR;
    double pL;

    switch(k)
    {
    case 10:
        for(int i = 0; i < pValue.size(); i++)
        {
            if(pSquare >= pStatic9[i])
            {
                pR = pValue[i];
                if(i != 0)
                    pL = pValue[i-1];
                else pL = 0;
                break;
            }
        }
        break;

    case 8:
        for(int i = 0; i < pValue.size(); i++)
        {
            if(pSquare >= pStatic7[i])
            {
                pR = pValue[i];
                if(i != 0)
                    pL = pValue[i-1];
                else pL = 0;
                break;
            }
        }
        break;

    case 14:
        for(int i = 0; i < pValue.size(); i++)
        {
            if(pSquare >= pStatic13[i])
            {
                pR = pValue[i];
                if(i != 0)
                    pL = pValue[i-1];
                else pL = 0;
                break;
            }
        }
        break;
    }

    cout << "p-value: " << pL << " <= p-value <= " << pR << endl;
    if (pR <= aValue)
        cout << "The hypothesis is rejected\n\n";
    else cout << "The hypothesis is accepted\n\n";
}

void pValueMMK(double pValue, double aValue = 0.05)
{
    if (pValue > aValue)
        cout << "The hypothesis is accepted";
    else cout << "The hypothesis is rejected";
}

int main()
{
    int n;
    cout << "Enter n: ";
    cin >> n;
    cout << "\n";

    //��� �������
    int k;
    cout << "Enter k: ";
    cin >> k;
    cout << "\n";

    char r;
    cout << "Choose distribution (r or e): ";
    cin >> r;
    cout << "\n";

    vector <double> sample;
    sample.resize(n);

    string sampleName;

    cout << "Enter file name for the sample: ";
    cin >> sampleName;
    cout << endl;

    //����� ����� � ��������
    ifstream Y(sampleName);
    if(!Y)
        {
            cout << "File not founded.\n\n";
            return 1;
        }

    //���������� ������� ��������
    for (int i = 0, j = 0; i < n; i++)
    {
        Y >> sample[i];
    }

    double xPSquare = pirsonSquare(sample, k, r);
    cout << "Pirson Square is " << xPSquare << "\n\n";

    pValueTab(xPSquare, k);

    double pValue = monteKarlo(xPSquare, n, k, r);
    cout << "MMK p-value is " << pValue << "\n";

    pValueMMK(pValue);

    cout << endl;

    return 0;
}
