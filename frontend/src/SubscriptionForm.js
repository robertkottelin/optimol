import React, { useState } from "react";
import { useStripe, useElements, CardElement } from "@stripe/react-stripe-js";
import axios from "axios";

const SubscriptionForm = ({ onSuccess }) => {
  const stripe = useStripe();
  const elements = useElements();
  const [email, setEmail] = useState(""); // State to store email input

  const handleSubmit = async (event) => {
    event.preventDefault();

    const cardElement = elements.getElement(CardElement);

    if (!email) {
      alert("Please enter your email.");
      return;
    }

    try {
      // Create a payment method
      const { error, paymentMethod } = await stripe.createPaymentMethod({
        type: "card",
        card: cardElement,
        billing_details: { email },
      });

      if (error) {
        console.error("Payment method error:", error);
        return;
      }

      // Send payment method ID and email to the backend
      const response = await axios.post("http://localhost:5000/subscribe", {
        email,
        paymentMethodId: paymentMethod.id,
      });

      if (response.data.success) {
        const { clientSecret } = response.data;

        // Confirm the subscription on the client side
        const { error: confirmError } = await stripe.confirmCardPayment(
          clientSecret
        );

        if (!confirmError) {
          console.log("Subscription successful");
          onSuccess(email); // Pass the email to parent component on success
        } else {
          console.error("Error confirming subscription:", confirmError);
        }
      } else {
        console.error("Subscription error:", response.data.error);
      }
    } catch (err) {
      console.error("Error in subscription process:", err);
    }
  };

  return (
    <div className="subscription-form">
      <h3>Subscribe for $10/month to get access to full compute power.</h3>
      <form onSubmit={handleSubmit}>
        <input
          type="email"
          value={email}
          onChange={(e) => setEmail(e.target.value)}
          placeholder="Enter your email"
          required
        />
        <div className="stripe-card-element">
          <CardElement
            options={{
              style: {
                base: {
                  color: "#fff",
                  fontSize: "16px",
                  "::placeholder": { color: "#888" },
                },
                invalid: { color: "#f44" },
              },
            }}
          />
        </div>
        <button type="submit" disabled={!stripe}>
          Subscribe
        </button>
      </form>
    </div>
  );
};

export default SubscriptionForm;
