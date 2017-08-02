class Body:
    """
    Class for the definition of a body with a given mass
    """
    def __init__(self, name, mass, position, velocity):
        self.name = name
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.acceleration = 0
        
    def position(self, position):
        """
        method for updating the velocity of the body
        """
        self.position = position
        
    def velocity(self, velocity):
        """
        Method for updating the position of the body
        """
        self.velocity = velocity

    def acceleration(self, acceleration):
        """
        Method for updating the acceleration of the
        body
        """
        self.acceleration = acceleration

    def printStatus(self):
        """
        Prints all the values for this body
        """
        print self.name
        print self.mass
        print self.position
        print self.velocity
        print self.acceleration
    
