// AuthContext.js
import React, { createContext, useState, useEffect } from 'react';
import axios from 'axios';

export const AuthContext = createContext({
  currentUser: null,
  isAuthenticated: false,
  isSubscribed: false,
  isLoading: true,
  login: () => {},
  register: () => {},
  logout: () => {}
});

export const AuthProvider = ({ children }) => {
  const [currentUser, setCurrentUser] = useState(null);
  const [isLoading, setIsLoading] = useState(true);
  const apiBaseUrl = "https://optimizemolecule.com";

  // Configure axios to include credentials with all requests
  axios.defaults.withCredentials = true;

  // Check if user is already logged in on mount
  useEffect(() => {
    const checkAuthStatus = async () => {
      try {
        const response = await axios.get(`${apiBaseUrl}/me`);
        setCurrentUser(response.data);
      } catch (error) {
        console.log("User not authenticated");
        setCurrentUser(null);
      } finally {
        setIsLoading(false);
      }
    };

    checkAuthStatus();
  }, []);

  // Login function
  const login = async (email, password) => {
    try {
      const response = await axios.post(`${apiBaseUrl}/login`, { email, password });
      setCurrentUser(response.data.user);
      return { success: true };
    } catch (error) {
      return { 
        success: false, 
        error: error.response?.data?.error || "Login failed" 
      };
    }
  };

  // Register function
  const register = async (email, password) => {
    try {
      const response = await axios.post(`${apiBaseUrl}/register`, { email, password });
      setCurrentUser(response.data.user);
      return { success: true };
    } catch (error) {
      return { 
        success: false, 
        error: error.response?.data?.error || "Registration failed" 
      };
    }
  };

  // Logout function
  const logout = async () => {
    try {
      await axios.post(`${apiBaseUrl}/logout`);
      setCurrentUser(null);
      return { success: true };
    } catch (error) {
      return { success: false, error: "Logout failed" };
    }
  };

  const value = {
    currentUser,
    isLoading,
    login,
    register,
    logout,
    isAuthenticated: !!currentUser,
    isSubscribed: currentUser?.isSubscribed || false
  };

  return <AuthContext.Provider value={value}>{children}</AuthContext.Provider>;
};