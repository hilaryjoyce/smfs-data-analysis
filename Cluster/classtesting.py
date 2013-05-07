class Top(object):
    
    def __init__(self, p=0.4):
        super(Top, self).__init__()
        self.p = p
        print "Top is here!"

class Above(Top):

    def __init__(self, s = 4):
        super(Above, self).__init__()
        self.s = s
        print "above all?"
