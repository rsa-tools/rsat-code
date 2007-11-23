

class MemoryCache(object):
    """
    simple memory cache
    key value association
    """
    def __init__(self, f):
        self.data = {}
        self.f = f
        #self.access = 0
        #self.success = 0
        
    def __call__(self, key):
        #self.access += 1
        if not self.data.has_key(key):
            self.data[key] = self.f(key)
        #else:
        #    self.success += 1
        return self.data[key]
        
    def stats(self):
        pass
        #return 'efficency: %.3f' % ( self.success / float(self.access))