
# coding: utf-8

# In[5]:


import hashlib

def geohash(latitude, longitude, datedow):
    # http://xkcd.com/426/
    h = hashlib.md5(datedow).hexdigest()
    p, q = [('%f' % float.fromhex('0.' + x)) for x in (h[:16], h[16:32])]
    print('%d%s %d%s' % (latitude, p[1:], longitude, q[1:]))


# In[6]:


geohash(47.0486395,-122.9229211,b'2018-02-20-24964.75')


# nope. that's underwater.
