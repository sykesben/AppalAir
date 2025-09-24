import pickle
import pyAesCrypt
import os

class MPLNetAccess():
    """
    Provides access to the MPLNet website.
    
    Attributes:
        username: A string representing the username.
        password: A string representing the password.
        
    Methods:
        createAccess: Creates a file containing the username and password.
        getAccess: Retrieves the username and password from the file.
        getUser: Returns the username.
        getPass: Returns the password.
        findDat: Finds the file containing the username and password.
    """

    def __init__(self):
        self.__username = ''
        self.__password = ''

    def createAccess(self, username, password):
        self.__username = username
        self.__password = password
        KeyAccess = "MPLNetAccess.key(AES@2023J)"
        data = MPLNetAccess()
        data.username = username
        data.password = password
        # pickle
        with open('mplnet.pkl', 'wb') as f:
            pickle.dump(data, f)
        # encrypt
        with open('mplnet.pkl', 'rb') as f:
            with open('mplnet.dat', 'wb') as g:
                pyAesCrypt.encryptStream(f, g, KeyAccess, 64*1024)
        # clean up unencrypted file
        os.remove('mplnet.pkl')

    def getAccess(self):
        """
        Retrieves the username and password from the file.
        
        Args:
            None
        Returns:
            True if successful, False otherwise.
        """
        KeyAccess = "MPLNetAccess.key(AES@2023J)"

        # find file
        file = self.findDat()
        if file == None:
            return False

        # decrypt
        # with open('mplnet.dat', 'rb') as f:
        with open(file, 'rb') as f:
            try:
                with open('mplnet.pkl', 'wb') as g:
                    pyAesCrypt.decryptStream(f, g, KeyAccess, 64*1024)
            except ValueError:
                print('Decryption failed.')
                os.remove('mplnet.pkl') 
        # unpickle
        with open('mplnet.pkl', 'rb') as f:
            data = pickle.load(f)
        # clean up unencrypted file
        os.remove('mplnet.pkl')
        self.__username = data.username
        self.__password = data.password

    def getUser(self):
        """
        Returns the username.
        
        Args:
            None
        Returns:
            A string representing the username.
        """

        return self.__username
    
    def getPass(self):
        """
        Returns the password.
        
        Args:
            None
        Returns:
            A string representing the password.
        """

        return self.__password

    def findDat(self):
        """
        Finds the encrypted file containing the username and password.
        
        Args:
            None
        Returns:
            A string representing the path to the file.
        """

        start = '..\\'
        file = 'mplnet.dat'
        for dir, dname, fname in sorted(os.walk(start)):
            for f in fname:
                if f == file:
                    return os.path.join(dir, f)
        return None
    
    