#include "provided.h"
#include <unordered_map>
#include "Trie.h"
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
using namespace std;

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    int m_minSearchLength;
    vector<Genome> genomes;
    Trie<pair<int,int>> trie;   //the pair is <name_index,position_index>
    //helper function of findGenomesWithThisDNA
    int findHelper(Genome gen,bool exactMatchOnly,int pos,const string& fragment) const
    {
        string previousFragment;
        gen.extract(pos, m_minSearchLength, previousFragment);
        //if the previousFragment are already a snip
        if (previousFragment!=fragment.substr(0,m_minSearchLength))
            exactMatchOnly=true;
        int i=m_minSearchLength;
        string thisChar;
        for (;i<fragment.length() && i+pos< gen.length();i++)
        {
            gen.extract(i+pos, 1, thisChar);
            //if the current character does not equal fragment’s corresponding character
            if (thisChar!=fragment.substr(i,1))
            {
                if (exactMatchOnly)
                    return i;
                exactMatchOnly=true;//give it a chance
            }
        }
        return i;
    }
    /*
    bool compareTwoGenomeMatch(const GenomeMatch& GM1, const GenomeMatch& GM2)
    {
        if (GM1.percentMatch>GM2.percentMatch)
            return true;
        if (GM1.percentMatch<GM2.percentMatch)
            return false;
        return (GM1.genomeName<GM2.genomeName);
    }
    bool compare(const GenomeMatch* GM1, const GenomeMatch* GM2)
    {
        return compareTwoGenomeMatch(*GM1, *GM2);
    }
     */
};

//set up
GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
    m_minSearchLength=minSearchLength;
}

//used to add a new genome to the library of genomes maintained by your GenomeMatcher object.
void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    genomes.push_back(genome);
    string fragment;
    for (int i=0;genome.extract(i, m_minSearchLength, fragment);i++)
    {
        //the pair is <name_index,position_index>
        trie.insert(fragment,make_pair(static_cast<int>(genomes.size()-1), i));
    }
}

//get minimum search length
int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minSearchLength;
}

//ued to find all genomes in the library that contain a specified DNA fragment (e.g., “GATTACA”), or potentially one or more of its SNiPs (e.g. “GCTTACA”, “GATTATA”), which are minimumLength or more bases long.
bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    std::vector<pair<int,int>> searchResult;
    unordered_map<std::string,pair<int,int>> tempMatches;
    bool found=false;
    /*
     The findGenomesWIthThisDNA() method must return false if
        1. fragment's length is less than minimumLength, or
        2. minimumLength is less than the minSearchLength value passed to the GnomeMatcher's constructor,
     */
    if (fragment.length()<minimumLength)
        return false;
    if (minimumLength<minimumSearchLength())
        return false;
    searchResult=trie.find(fragment.substr(0,m_minSearchLength), exactMatchOnly);
    int totalLength;
    for (size_t k=0;k<searchResult.size();k++)
    {
        //use findHelper to get the total length with current search result’s genome’s name index
        totalLength=findHelper(genomes[searchResult[k].first],exactMatchOnly,searchResult[k].second,fragment);
        if (totalLength<minimumLength)
            continue;
        //check if such search result’s genome’s name is already in the tempMatches
        auto search=tempMatches.find(genomes[searchResult[k].first].name());
        //if so, update the tempMatches with the longer total length
        if (search != tempMatches.end())
        {
            if (search->second.first<totalLength)
            {
                search->second.first=totalLength;
                search->second.second=searchResult[k].second;
            }
            continue;
        }
        //if not, insert the new search result’s genome name and total length and position to the tempMatches
        tempMatches.insert(make_pair(genomes[searchResult[k].first].name(), make_pair(totalLength,searchResult[k].second)));
        found=true;
    }
    if (!found)
        return false;
    //push back all the items in the unordered_map tempMatches to vector matches
    unordered_map<std::string,pair<int,int>>:: iterator itr;
    for (itr = tempMatches.begin(); itr != tempMatches.end(); itr++)
    {
        DNAMatch thisDNA;
        thisDNA.genomeName=itr->first;
        thisDNA.length=itr->second.first;
        thisDNA.position=itr->second.second;
        matches.push_back(thisDNA);
    }
    return true;
}
/*
bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    if (fragmentMatchLength<minimumSearchLength())
        return false;
    bool found=false;
    int S=query.length()/fragmentMatchLength;
    //int count=0;
    vector<GenomeMatch> totalMatches;
    for (int i=0;i<S;i++)
    {
        string partSequence;
        query.extract(i*fragmentMatchLength, fragmentMatchLength, partSequence);
        vector<DNAMatch> currentMatches;
        findGenomesWithThisDNA(partSequence, fragmentMatchLength, exactMatchOnly, currentMatches);
        //if(findGenomesWithThisDNA(partSequence, fragmentMatchLength, exactMatchOnly, currentMatches)) count+=currentMatches.size();
        for (size_t k=0;k<currentMatches.size();k++)
        {
            found=true;
            bool addNew=true;
            for (size_t j=0;j<totalMatches.size();j++)
            {
                if (totalMatches[j].genomeName==currentMatches[k].genomeName)
                {
                    totalMatches[j].percentMatch+=100.00/S;
                    addNew=false;
                    break;
                }
            }
            if (!addNew)
                continue;
            GenomeMatch thisGenomeMatch;
            thisGenomeMatch.genomeName=currentMatches[k].genomeName;
            thisGenomeMatch.percentMatch=100.00/S;
            totalMatches.push_back(thisGenomeMatch);
        }
    }
    for (size_t j=0;j<totalMatches.size();j++)
    {
        if (totalMatches[j].percentMatch>=matchPercentThreshold)
            results.push_back(totalMatches[j]);
    }
    sort(results.begin(),
         results.end(),
         [](const GenomeMatch& lhs, const GenomeMatch& rhs)
    {
        if (lhs.percentMatch>rhs.percentMatch)
            return true;
        if (lhs.percentMatch<rhs.percentMatch)
            return false;
        return (lhs.genomeName<rhs.genomeName);
    });
    return found;
}
*/


//The findRelatedGenomes() method compares a passed-in query genome for a new organism against all genomes currently held in a GenomeMatcher object’s library and passes back a vector of all genomes that contain more than matchPercentThreshold of the base sequences of length fragmentMatchLength from the query genome.
bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    if (fragmentMatchLength<minimumSearchLength())
        return false;
    bool found=false;
    int S=query.length()/fragmentMatchLength;
    unordered_map<std::string,double> totalMatches;
    for (int i=0;i<S;i++)
    {
        string partSequence;
        //Extract that sequence from the queried genome.
        query.extract(i*fragmentMatchLength, fragmentMatchLength, partSequence);
        vector<DNAMatch> currentMatches;
        //Search for the extracted sequence across all genomes in the library (using findGenomesWithThisDNA())
        findGenomesWithThisDNA(partSequence, fragmentMatchLength, exactMatchOnly, currentMatches);
        //If a match is found in one or more genomes in the library, then for each such genome, increase the count of matches found thus far for it.
        for (size_t k=0;k<currentMatches.size();k++)
        {
            found=true;
            bool addNew=true;
            for (size_t j=0;j<totalMatches.size();j++)
            {
                auto search = totalMatches.find(currentMatches[k].genomeName);
                if (search != totalMatches.end())
                {
                    search->second+=100.00/S;
                    addNew=false;
                    break;
                }
            }
            if (!addNew)
                continue;
            totalMatches.insert(make_pair(currentMatches[k].genomeName, 100.00/S));
        }
    }
    if (!found)
        return false;
    //push back all totalmatches to results
    unordered_map<std::string, double>:: iterator itr;
    for (itr = totalMatches.begin(); itr != totalMatches.end(); itr++)
    {
        if (itr->second>=matchPercentThreshold)
        {
            GenomeMatch thisGenomeMatch;
            thisGenomeMatch.genomeName=itr->first;
            thisGenomeMatch.percentMatch=itr->second;
            results.push_back(thisGenomeMatch);
        }
    }
    //ordered in descending order by the match proportion p, and breaking ties by the genome name in ascending alphabetical order.
    sort(results.begin(),
         results.end(),
         [](const GenomeMatch& lhs, const GenomeMatch& rhs)
         {
             if (lhs.percentMatch>rhs.percentMatch)
                 return true;
             if (lhs.percentMatch<rhs.percentMatch)
                 return false;
             return (lhs.genomeName<rhs.genomeName);
         });
    return found;
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
