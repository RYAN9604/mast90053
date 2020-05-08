#!/usr/bin/env python
# coding: utf-8

# In[2]:


k = var("k")


# In[3]:


t = factorial(k)


# In[4]:


t


# In[5]:


t.simplify_full()


# In[7]:


f = t(k+1)/t(k)


# In[9]:


fs = f.simplify_full()


# In[10]:


fs


# In[11]:


type(fs)


# In[12]:


fs.is_polynomial(k)


# In[13]:


n = var("n")
t = binomial(n, k)


# In[14]:


t


# In[15]:


f = t(k+1)/t(k)


# In[16]:


f


# In[17]:


fs = f.simplify_full()


# In[18]:


fs


# In[19]:


num, denom = fs.numerator_denominator()


# In[20]:


num


# In[21]:


denom


# In[22]:


num.is_polynomial(k)


# In[23]:


denom.is_polynomial(k)


# In[24]:


t = (-1)^k * x^(2*k+1) / factorial(2*k+1)


# In[25]:


t


# In[26]:


show(t)


# In[28]:


f = t(k+1)/t(k)


# In[29]:


f


# In[30]:


show(f)


# In[31]:


fs = f.simplify_full()


# In[32]:


fs


# In[33]:


show(fs)


# In[34]:


num, denom = fs.numerator_denominator()


# In[35]:


num.is_polynomial(k)


# In[40]:


def is_hypergeometric_term(t, k):
    ratio = t.subs({k: k+1}) / t.subs({k: k})
    ratio = ratio.simplify_full()
    num, denom = ratio.numerator_denominator()
    return (num.is_polynomial(k) and denom.is_polynomial(k))


# In[41]:


t


# In[42]:


is_hypergeometric_term(t, k)


# In[43]:


t = factorial(k)


# In[44]:


is_hypergeometric_term(t, k)


# In[45]:


is_hypergeometric_term(log(k^2-tan(k)), k)


# In[ ]:




