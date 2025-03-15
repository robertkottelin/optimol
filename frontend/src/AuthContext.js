import React, { createContext, useState, useEffect } from 'react';
import axios from 'axios';

export const AuthContext = createContext({
  currentUser: null,
  isAuthenticated: false,
  isSubscribed: false,
  isLoading: true,
  error: null,
  login: () => Promise.resolve({ success: false }),
  register: () => Promise.resolve({ success: false }),
  logout: () => Promise.resolve({ success: false })
});

export const AuthProvider = ({ children }) => {
  const [state, setState] = useState({
    currentUser: null,
    isAuthenticated: false,
    isSubscribed: false,
    isLoading: true,
    error: null
  });
  
  const apiBaseUrl = "https://optimizemolecule.com";

  // Configure axios defaults
  axios.defaults.withCredentials = true;

  // Handle authentication status check
  useEffect(() => {
    let mounted = true;
    
    const checkAuth = async () => {
      try {
        const response = await axios.get(`${apiBaseUrl}/me`);
        if (!mounted) return;
        
        setState(prev => ({
          ...prev,
          currentUser: response.data,
          isAuthenticated: true,
          isSubscribed: response.data.isSubscribed || false,
          isLoading: false
        }));
      } catch (error) {
        if (!mounted) return;
        
        // Expected error for unauthenticated users - not a failure case
        if (error.response?.status === 401) {
          setState(prev => ({
            ...prev,
            currentUser: null,
            isAuthenticated: false,
            isSubscribed: false,
            isLoading: false,
            error: null
          }));
          return;
        }
        
        // Unexpected error
        console.error("Authentication check failed:", error);
        setState(prev => ({
          ...prev,
          currentUser: null,
          isAuthenticated: false,
          isSubscribed: false,
          isLoading: false,
          error: error.message || "Authentication check failed"
        }));
      }
    };

    checkAuth();
    return () => { mounted = false; };
  }, []);

  // Login function
  const login = async (email, password) => {
    try {
      setState(prev => ({ ...prev, isLoading: true }));
      const response = await axios.post(`${apiBaseUrl}/login`, { email, password });
      
      setState(prev => ({
        ...prev,
        currentUser: response.data.user,
        isAuthenticated: true,
        isSubscribed: response.data.user.isSubscribed || false,
        isLoading: false,
        error: null
      }));
      
      return { success: true };
    } catch (error) {
      setState(prev => ({ 
        ...prev, 
        isLoading: false,
        error: error.response?.data?.error || "Login failed"
      }));
      
      return { 
        success: false, 
        error: error.response?.data?.error || "Login failed" 
      };
    }
  };

  // Register function
  const register = async (email, password) => {
    try {
      setState(prev => ({ ...prev, isLoading: true }));
      const response = await axios.post(`${apiBaseUrl}/register`, { email, password });
      
      setState(prev => ({
        ...prev,
        currentUser: response.data.user,
        isAuthenticated: true,
        isSubscribed: response.data.user.isSubscribed || false,
        isLoading: false,
        error: null
      }));
      
      return { success: true };
    } catch (error) {
      setState(prev => ({ 
        ...prev, 
        isLoading: false,
        error: error.response?.data?.error || "Registration failed"
      }));
      
      return { 
        success: false, 
        error: error.response?.data?.error || "Registration failed" 
      };
    }
  };

  // Logout function
  const logout = async () => {
    try {
      setState(prev => ({ ...prev, isLoading: true }));
      await axios.post(`${apiBaseUrl}/logout`);
      
      setState({
        currentUser: null,
        isAuthenticated: false,
        isSubscribed: false,
        isLoading: false,
        error: null
      });
      
      return { success: true };
    } catch (error) {
      setState(prev => ({ 
        ...prev, 
        isLoading: false,
        error: error.response?.data?.error || "Logout failed"
      }));
      
      return { 
        success: false, 
        error: error.response?.data?.error || "Logout failed" 
      };
    }
  };

  return (
    <AuthContext.Provider value={{
      ...state,
      login,
      register,
      logout
    }}>
      {children}
    </AuthContext.Provider>
  );
};