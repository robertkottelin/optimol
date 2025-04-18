import React, { createContext, useState, useEffect } from 'react';
import axios from 'axios';

// Enhanced JWT validation utility function
const isValidJWT = (token) => {
  if (!token) return false;
  
  try {
    // Basic format validation
    const tokenParts = token.split('.');
    if (tokenParts.length !== 3) return false;
    
    // Decode payload (won't validate signature - server handles that)
    const base64Url = tokenParts[1];
    const base64 = base64Url.replace(/-/g, '+').replace(/_/g, '/');
    const payload = JSON.parse(window.atob(base64));
    
    // We don't validate sub type here - we'll let server rejection handle it
    return true;
  } catch (e) {
    console.error("JWT validation error:", e);
    return false;
  }
};

export const AuthContext = createContext({
  currentUser: null,
  isAuthenticated: false,
  isSubscribed: false,
  isLoading: true,
  error: null,
  token: null,
  login: () => Promise.resolve({ success: false }),
  register: () => Promise.resolve({ success: false }),
  registerAndSubscribe: () => Promise.resolve({ success: false }),
  logout: () => Promise.resolve({ success: false })
});

export const AuthProvider = ({ children }) => {
  const [state, setState] = useState({
    currentUser: null,
    isAuthenticated: false,
    isSubscribed: false,
    isLoading: true,
    error: null,
    token: localStorage.getItem('access_token') || null
  });
  
  const apiBaseUrl = "/api";

  // Configure axios defaults
  useEffect(() => {
    // Configure axios defaults
    axios.defaults.headers.common['X-Requested-With'] = 'XMLHttpRequest';
    
    // Add token to headers via interceptor
    const interceptor = axios.interceptors.request.use(
      config => {
        if (state.token) {
          config.headers['Authorization'] = `Bearer ${state.token}`;
        }
        return config;
      },
      error => Promise.reject(error)
    );
    
    // Cleanup interceptor on unmount
    return () => axios.interceptors.request.eject(interceptor);
  }, [state.token]);

  // Handle authentication status check
  useEffect(() => {
    let mounted = true;
    
    const checkAuth = async () => {
      // Validate token format before attempting to use it
      if (!state.token || !isValidJWT(state.token)) {
        localStorage.removeItem('access_token');
        setState(prev => ({
          ...prev,
          token: null,
          isAuthenticated: false,
          isLoading: false
        }));
        return;
      }
      
      try {
        const response = await axios({
          method: 'get',
          url: `${apiBaseUrl}/me`,
          headers: {
            'Accept': 'application/json',
            'Content-Type': 'application/json'
          }
        });
        
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
        
        // Specific handling for "Subject must be a string" error
        if (error.response?.data?.msg === "Subject must be a string") {
          localStorage.removeItem('access_token');
          setState(prev => ({
            ...prev,
            token: null,
            currentUser: null,
            isAuthenticated: false,
            isSubscribed: false,
            isLoading: false,
            error: null
          }));
          return;
        }
        
        // Handle both 401 (Unauthorized) and 422 (Unprocessable Entity) - both indicate token issues
        if (error.response?.status === 401 || error.response?.status === 422) {
          // Clear invalid token
          localStorage.removeItem('access_token');
          
          setState(prev => ({
            ...prev,
            token: null,
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
  }, [state.token, apiBaseUrl]);

  // Login function
  const login = async (email, password) => {
    try {
      setState(prev => ({ ...prev, isLoading: true }));
      const response = await axios({
        method: 'post',
        url: `${apiBaseUrl}/login`,
        data: { email, password },
        headers: {
          'Content-Type': 'application/json',
          'Accept': 'application/json'
        }
      });
      
      // Store token in localStorage
      if (response.data.token) {
        localStorage.setItem('access_token', response.data.token);
      }
      
      setState(prev => ({
        ...prev,
        token: response.data.token,
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
      const response = await axios({
        method: 'post',
        url: `${apiBaseUrl}/register`,
        data: { email, password },
        headers: {
          'Content-Type': 'application/json',
          'Accept': 'application/json'
        }
      });
      
      // Store token in localStorage
      if (response.data.token) {
        localStorage.setItem('access_token', response.data.token);
      }
      
      setState(prev => ({
        ...prev,
        token: response.data.token,
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

  // New function: Register and subscribe in one step
  const registerAndSubscribe = async (email, password, paymentMethodId) => {
    try {
      setState(prev => ({ ...prev, isLoading: true }));
      
      // Step 1: Register the user
      const registerResponse = await axios({
        method: 'post',
        url: `${apiBaseUrl}/register`,
        data: { email, password },
        headers: {
          'Content-Type': 'application/json',
          'Accept': 'application/json'
        }
      });
      
      if (!registerResponse.data.token) {
        throw new Error("Registration successful but no token received");
      }
      
      // Store token and update auth state
      const token = registerResponse.data.token;
      localStorage.setItem('access_token', token);
      
      // Step 2: Subscribe with the new account
      const subscribeResponse = await axios.post(
        `${apiBaseUrl}/subscribe`, 
        { paymentMethodId },
        { 
          headers: {
            'Content-Type': 'application/json',
            'Authorization': `Bearer ${token}`
          }
        }
      );
      
      if (!subscribeResponse.data.success) {
        throw new Error(subscribeResponse.data.error || "Subscription failed");
      }
      
      // Update auth state with the registered and subscribed user
      setState(prev => ({
        ...prev,
        token: token,
        currentUser: registerResponse.data.user,
        isAuthenticated: true,
        isSubscribed: true,
        isLoading: false,
        error: null
      }));
      
      return { 
        success: true, 
        clientSecret: subscribeResponse.data.clientSecret 
      };
    } catch (error) {
      setState(prev => ({ 
        ...prev, 
        isLoading: false,
        error: error.response?.data?.error || error.message || "Registration and subscription failed"
      }));
      
      return { 
        success: false, 
        error: error.response?.data?.error || error.message || "Registration and subscription failed"
      };
    }
  };

  // Logout function
  const logout = async () => {
    try {
      setState(prev => ({ ...prev, isLoading: true }));
      await axios.post(`${apiBaseUrl}/logout`);
      
      // Remove token from localStorage
      localStorage.removeItem('access_token');
      
      setState({
        currentUser: null,
        isAuthenticated: false,
        isSubscribed: false,
        isLoading: false,
        error: null,
        token: null
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
      registerAndSubscribe,
      logout
    }}>
      {children}
    </AuthContext.Provider>
  );
};