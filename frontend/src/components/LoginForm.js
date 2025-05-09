import React, { useState, useContext } from "react";
import { AuthContext } from "../AuthContext";
import { Icons } from "./Icons";

const LoginForm = ({ toggleForm, onAuthSuccess }) => { // Add onAuthSuccess parameter
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState("");

  const { login } = useContext(AuthContext);

  const handleSubmit = async (e) => {
    e.preventDefault();
    setError("");
    
    if (!email || !password) {
      setError("Please enter both email and password");
      return;
    }

    setIsLoading(true);

    try {
      const result = await login(email, password);
      if (!result.success) {
        setError(result.error);
      } else {
        // Call onAuthSuccess if authentication was successful
        if (onAuthSuccess) {
          onAuthSuccess();
        }
      }
    } catch (err) {
      setError("An unexpected error occurred. Please try again.");
      console.error("Login error:", err);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div className="login-form" style={{ padding: "30px" }}>
      <h3 style={{ 
        fontSize: "24px", 
        fontWeight: "700", 
        marginBottom: "20px", 
        display: "flex", 
        alignItems: "center", 
        gap: "10px" 
      }}>
        <Icons.verified /> Sign In
      </h3>

      {error && (
        <div style={{ 
          color: "#f43f5e", 
          marginBottom: "15px", 
          padding: "10px", 
          borderRadius: "6px",
          backgroundColor: "rgba(244, 63, 94, 0.1)",
          border: "1px solid rgba(244, 63, 94, 0.2)"
        }}>
          {error}
        </div>
      )}

      <form onSubmit={handleSubmit}>
        <div style={{ marginBottom: "15px" }}>
          <label 
            htmlFor="email" 
            style={{ 
              display: "block", 
              marginBottom: "5px", 
              fontSize: "14px", 
              fontWeight: "500", 
              color: "#94a3b8" 
            }}
          >
            Email Address
          </label>
          <div style={{ position: "relative" }}>
            <span style={{ 
              position: "absolute", 
              left: "12px", 
              top: "50%", 
              transform: "translateY(-50%)", 
              color: "#94a3b8" 
            }}>
              <Icons.info />
            </span>
            <input
              id="email"
              type="email"
              value={email}
              onChange={(e) => setEmail(e.target.value)}
              style={{ 
                width: "100%", 
                padding: "12px 12px 12px 36px", 
                backgroundColor: "rgba(15, 23, 42, 0.7)", 
                border: "1px solid rgba(255, 255, 255, 0.1)", 
                borderRadius: "6px", 
                color: "#f0f4f8", 
                fontSize: "14px" 
              }}
              placeholder="your@email.com"
              required
            />
          </div>
        </div>

        <div style={{ marginBottom: "20px" }}>
          <label 
            htmlFor="password" 
            style={{ 
              display: "block", 
              marginBottom: "5px", 
              fontSize: "14px", 
              fontWeight: "500", 
              color: "#94a3b8" 
            }}
          >
            Password
          </label>
          <div style={{ position: "relative" }}>
            <span style={{ 
              position: "absolute", 
              left: "12px", 
              top: "50%", 
              transform: "translateY(-50%)", 
              color: "#94a3b8" 
            }}>
              <Icons.lock />
            </span>
            <input
              id="password"
              type="password"
              value={password}
              onChange={(e) => setPassword(e.target.value)}
              style={{ 
                width: "100%", 
                padding: "12px 12px 12px 36px", 
                backgroundColor: "rgba(15, 23, 42, 0.7)", 
                border: "1px solid rgba(255, 255, 255, 0.1)", 
                borderRadius: "6px", 
                color: "#f0f4f8", 
                fontSize: "14px" 
              }}
              placeholder="••••••••"
              required
            />
          </div>
        </div>

        <button
          type="submit"
          disabled={isLoading}
          className="float"
          style={{
            width: "100%",
            padding: "12px",
            borderRadius: "8px",
            backgroundColor: "rgba(56, 189, 248, 0.9)",
            color: "white",
            fontSize: "16px",
            fontWeight: "600",
            border: "none",
            cursor: isLoading ? "not-allowed" : "pointer",
            opacity: isLoading ? "0.7" : "1",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            gap: "8px",
            marginBottom: "20px"
          }}
        >
          {isLoading ? (
            <>
              <span className="spinner"></span>
              Signing In...
            </>
          ) : (
            "Sign In"
          )}
        </button>

        <div style={{ textAlign: "center", fontSize: "14px", color: "#94a3b8" }}>
          Don't have an account?{" "}
          <button
            type="button"
            onClick={toggleForm}
            style={{
              background: "none",
              border: "none",
              color: "#38bdf8",
              cursor: "pointer",
              fontWeight: "500",
              padding: "0",
              textDecoration: "underline"
            }}
          >
            Register now
          </button>
        </div>
      </form>
    </div>
  );
};

export default LoginForm;