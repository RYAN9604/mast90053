#!/usr/bin/env python
# coding: utf-8

# In[1]:


k = var("k")
t = (-1)^k * x^(2*k+1) / factorial(2*k+1)


# In[2]:


show(t)


# In[4]:


f = t(k+1)/t(k)


# In[7]:


fs = f.simplify_full()


# In[6]:


show(_)


# In[8]:


fs


# In[9]:


fs.factor()


# In[10]:


show(_)


# In[11]:


# this is 0F1[-,3/2;(-x^2/4)]


# In[21]:


t = 1 / ((2*k+1) * gamma(2*k+4))


# In[22]:


t


# In[23]:


show(t)


# In[24]:


t(1)


# In[25]:


t(0)


# In[26]:


t(-1)


# In[27]:


t(-2)


# In[28]:


t(-3)


# In[29]:


u = t.subs({k: k-1})


# In[30]:


u


# In[31]:


show(u)


# In[33]:


f = u(k+1) / u(k)


# In[35]:


fs = f.simplify_full()


# In[36]:


fs.factor()


# In[37]:


show(_)


# In[38]:


# this is  - 1F2[-1/2, 3/2 1/2; 1/4]


# In[39]:


u(0)


# In[ ]:




