#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
    string m_name;
    string m_sequence;
};

//set up
GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
    m_name=nm;
    m_sequence=sequence;
}

//load to genomes from files
bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes)
{
    string line;
    string name;
    string sequence;
    //get the first line
    getline(genomeSource,line);
    if (line.empty())
        return false;
    if (line[0]!='>')
        return false;
    if (line.substr(1).empty())
        return false;
    name=line.substr(1);
    //while not end of the file
    while(getline(genomeSource,line))
    {
        //can't be empty line in the file
        if (line.empty())
            return false;
        //a name?
        if (line[0]=='>')
        {
            if (sequence.empty())
                return false;
            genomes.push_back(Genome(name,sequence));
            name=line.substr(1);
            sequence.clear();
        }
        //a dna sequence?
        else
        {
            //up to 80
            if (line.size()>80)
                return false;
            for (int i=0;i<line.size();i++)
            {
                if (line[i]!='A' && line[i]!='a' && line[i]!='G' && line[i]!='g' && line[i]!='T' && line[i]!='t' && line[i]!='C' && line[i]!='c' && line[i]!='N' && line[i]!='n')
                    return false;
            }
            sequence+=line;
        }
    }
    //the last line should also not be empty
    if (sequence.empty())
        return false;
    genomes.push_back(Genome(name,sequence));
    return true;
}


//return the length of sequence
int GenomeImpl::length() const
{
    return static_cast<int>(m_sequence.size());
}

//return the name
string GenomeImpl::name() const
{
    return m_name;
}


//The extract() method must set fragment to a copy of a portion of the Genome's DNA sequence: the substring length characters long starting at position (where the first character of the sequence is at position 0).
//The extract() method must return true if it successfully extracts a string of the specified length, and false otherwise
bool GenomeImpl::extract(int position, int length, string& fragment) const
{
    if ((position+length)>this->length())
        return false;
    fragment=m_sequence.substr(position,length);
    return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes)
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
