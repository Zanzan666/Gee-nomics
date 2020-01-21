#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>

template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie();
    void reset();
    void insert(const std::string& key, const ValueType& value);
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;
    
    // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
    struct TreeNode
    {
        std::vector <TreeNode*> m_childrenPtr;
        std::vector <char> m_childrenLabel;
        std::vector<ValueType> m_value;
    };
    TreeNode* m_root;
    //this is the helper function of insert
    TreeNode* insertHelper(std::string key,TreeNode* p)
    {
        //reached the end
        if (key=="")
            return p;
        //for all child pointers and the corresponding labels
        for (size_t k=0;k<p->m_childrenPtr.size();k++)
        {
            if (p->m_childrenLabel[k]==key[0])
                //move on to the next char in key with this childrenPtr
                return insertHelper(key.substr(1),p->m_childrenPtr[k]);
        }
        //no such label exists, create new one
        TreeNode* thisChildren=new TreeNode;
        p->m_childrenPtr.push_back(thisChildren);
        p->m_childrenLabel.push_back(key[0]);
        //move on to the next char in key with this childrenPtr
        return insertHelper(key.substr(1),p->m_childrenPtr.back());
    }
    /*
    std::vector<ValueType> findHelper(std::string key, bool canBeWrong, TreeNode* p) const
    {
        if (p==nullptr)
        {
            std::vector<ValueType> empty;
            return empty;
        }
        if (key=="")
            return p->m_value;
        std::vector<ValueType> totalValue;
        for (size_t k=0;k<p->m_childrenPtr.size();k++)
        {
            std::vector<ValueType> thisValue;
            if (p->m_childrenLabel[k]==key[0])
                thisValue=findHelper(key.substr(1), canBeWrong, p->m_childrenPtr[k]);
            else if (canBeWrong && p!=m_root)
                thisValue=findHelper(key.substr(1), false, p->m_childrenPtr[k]);
            for (size_t k=0;k<thisValue.size();k++)
                totalValue.push_back(thisValue[k]);
            
        }
        return totalValue;
    }
    */
    void findHelper(std::string key, bool canBeWrong, TreeNode* p,std::vector<ValueType>& result) const
    {
        //reaches the end of the tree
        if (p==nullptr)
            return;
        //reaches the end of the key, therefore we can put in the values at the current position to the result
        if (key.empty())
        {
            for (size_t k=0;k<p->m_value.size();k++)
                result.push_back(p->m_value[k]);
            return;
        }
        //for all child pointers and the correpsonding labels
        for (size_t k=0;k<p->m_childrenPtr.size();k++)
        {
            //if matches
            if (p->m_childrenLabel[k]==key[0])
                //move on to the next
                findHelper(key.substr(1), canBeWrong, p->m_childrenPtr[k],result);
            //if not matches and can be wrong and is not at the beginning
            else if (canBeWrong && p!=m_root)
                //move on to the next with can't be wrong
                findHelper(key.substr(1), false, p->m_childrenPtr[k],result);
            
        }
    }
    
    //recursivly free all dynamically allocated pointers in the trie
    void cleaner(TreeNode* p)
    {
        if (p==nullptr)
            return;
        for (size_t k=0;k<p->m_childrenPtr.size();k++)
            cleaner( p->m_childrenPtr[k]);
        delete p;
    }
};

//set up by create a root
template<typename ValueType>
Trie<ValueType>::Trie()
{
    m_root=new TreeNode;
}

//destruct by calling the cleaner
template<typename ValueType>
Trie<ValueType>::~Trie()
{
    cleaner(m_root);
}

//first clean and then allocate a new root
template<typename ValueType>
void Trie<ValueType>::reset()
{
    cleaner(m_root);
    m_root=new TreeNode;
}

//call the insert helper function and the returned pointer is where the value should store
template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value)
{
    TreeNode* StoreValueHere;
    StoreValueHere=insertHelper(key, m_root);
    StoreValueHere->m_value.push_back(value);
}

//call the find helper function
template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
    std::vector<ValueType> result;
    findHelper(key,!exactMatchOnly,m_root,result);
    return result;
    //return findHelper(key,!exactMatchOnly,m_root);
}



#endif // TRIE_INCLUDED

